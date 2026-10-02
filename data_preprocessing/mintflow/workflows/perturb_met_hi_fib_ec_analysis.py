# Scientific workflow derived from Revision/Nan annotation/Data/RMD files/MintFlow/00_scripts_and_reproducibility/perturb_met_hi_fib_ec_analysis.py
import os
from pathlib import Path
def project_path(relative):
    return str(Path(os.environ.get("METABOLIC_INPUT_DIR", "inputs")) / relative)
#!/usr/bin/env python

import argparse
import html
import json
import os
import pickle
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import squidpy as sq
import torch
from scipy import sparse
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

import mintflow


EC_LABELS = ["Met_hi_EC", "Met_int_EC"]
SCORE_GENESETS = {
    "glycolysis_score": ["SLC2A1", "GLUT1", "HK1", "PFKL", "PFKM", "PKM", "PKM2", "LDHA", "GAPDH"],
    "tca_score": ["CS", "OGDH", "ATP5F1A", "ATP5A", "SDHA"],
    "hif1a_score": ["HIF1A", "HIF1a"],
    "nox4_score": ["NOX4"],
}


def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)


def make_spatial_graph(adata, n_neighs=5):
    adata.uns = {}
    sq.gr.spatial_neighbors(
        adata=adata,
        spatial_key="spatial",
        library_key=None,
        set_diag=False,
        delaunay=False,
        n_neighs=n_neighs,
    )


def avg_generated_xmic(result):
    return np.stack(
        [r["MintFLow_Generated_Xmic"] for r in result["list_generated_realisations_ie_expressions"]],
        axis=0,
    ).mean(axis=0)


def as_log1p_nonnegative(X):
    if sparse.issparse(X):
        X = X.copy()
        X.data[X.data < 0] = 0
        X.data = np.log1p(X.data)
        return X
    return np.log1p(np.clip(np.asarray(X), 0, None))


def train_ec_classifier(adata):
    ec_mask = adata.obs["metEClabel"].isin(EC_LABELS).to_numpy()
    X_train = as_log1p_nonnegative(adata.X[ec_mask])
    y_train = adata.obs.loc[ec_mask, "metEClabel"].astype(str).to_numpy()
    clf = make_pipeline(
        StandardScaler(with_mean=False),
        LogisticRegression(max_iter=5000, class_weight="balanced", solver="liblinear"),
    )
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_train)
    report = classification_report(y_train, y_pred, output_dict=True, zero_division=0)
    return clf, report


def train_ec_classifier_from_generated(adata, Xmic):
    ec_mask = adata.obs["metEClabel"].isin(EC_LABELS).to_numpy()
    X_train = as_log1p_nonnegative(Xmic[ec_mask, :])
    y_train = adata.obs.loc[ec_mask, "metEClabel"].astype(str).to_numpy()
    clf = make_pipeline(
        StandardScaler(with_mean=False),
        LogisticRegression(max_iter=5000, class_weight="balanced", solver="liblinear"),
    )
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_train)
    report = classification_report(y_train, y_pred, output_dict=True, zero_division=0)
    return clf, report


def find_gene_indices(var, genes):
    upper_to_idx = {}
    for i, (idx, row) in enumerate(var.iterrows()):
        candidates = [idx, row.get("gene_id", ""), row.get("gene_symbol", "")]
        for candidate in candidates:
            if pd.notna(candidate):
                upper_to_idx[str(candidate).upper()] = i
    found = []
    for gene in genes:
        idx = upper_to_idx.get(gene.upper())
        if idx is not None and idx not in found:
            found.append(idx)
    return found


def compute_scores(X, var):
    X = as_log1p_nonnegative(X)
    score_df = pd.DataFrame(index=np.arange(X.shape[0]))
    used_genes = {}
    for score_name, genes in SCORE_GENESETS.items():
        idx = find_gene_indices(var, genes)
        used_genes[score_name] = [str(var.iloc[i].get("gene_symbol", var.index[i])) for i in idx]
        if len(idx) == 0:
            score_df[score_name] = np.nan
            continue
        sub = X[:, idx]
        vals = np.asarray(sub.mean(axis=1)).ravel()
        score_df[score_name] = vals
    return score_df, used_genes


def mcc_changed_by_cell_id(result_base, adata_base, result_scenario, adata_scenario):
    base_ids = adata_base.obs["cell_id"].astype(str).tolist()
    scenario_ids = adata_scenario.obs["cell_id"].astype(str).tolist()
    base_index = {cid: i for i, cid in enumerate(base_ids)}
    scenario_index = {cid: i for i, cid in enumerate(scenario_ids)}
    changed = {}
    common = set(base_index).intersection(scenario_index)
    for cid in common:
        changed[cid] = not np.allclose(
            result_base["np_MCC"][base_index[cid]],
            result_scenario["np_MCC"][scenario_index[cid]],
        )
    return changed


def analyze_scenario(name, adata_scenario, Xmic, clf, adata_base, result_base=None, result_scenario=None):
    ec_mask = adata_scenario.obs["metEClabel"].isin(EC_LABELS).to_numpy()
    obs = adata_scenario.obs.loc[ec_mask, ["cell_id", "sample_id", "metEClabel", "combined_met_label"]].copy()
    obs = obs.reset_index(drop=True)
    obs["scenario"] = name

    X_ec = Xmic[ec_mask, :]
    X_ec_log = as_log1p_nonnegative(X_ec)
    obs["predicted_ec_state"] = clf.predict(X_ec_log)
    classes = list(clf.named_steps["logisticregression"].classes_)
    prob = clf.predict_proba(X_ec_log)
    obs["prob_Met_hi_EC"] = prob[:, classes.index("Met_hi_EC")]

    scores, used_genes = compute_scores(X_ec, adata_scenario.var)
    obs = pd.concat([obs, scores], axis=1)

    obs["mcc_changed_vs_baseline"] = False
    if result_base is not None and result_scenario is not None:
        changed = mcc_changed_by_cell_id(result_base, adata_base, result_scenario, adata_scenario)
        obs["mcc_changed_vs_baseline"] = obs["cell_id"].astype(str).map(changed).fillna(False).astype(bool)

    return obs, used_genes


def summarise_cells(cell_df):
    score_cols = list(SCORE_GENESETS.keys())
    rows = []
    for changed_only, df in [("all_ecs", cell_df), ("mcc_changed_ecs", cell_df[cell_df["mcc_changed_vs_baseline"]])]:
        if df.empty:
            continue
        g = df.groupby(["scenario", "metEClabel"], observed=True)
        summary = g.agg(
            n_ec=("cell_id", "size"),
            n_pred_Met_hi_EC=("predicted_ec_state", lambda s: int((s == "Met_hi_EC").sum())),
            freq_pred_Met_hi_EC=("predicted_ec_state", lambda s: float((s == "Met_hi_EC").mean())),
            mean_prob_Met_hi_EC=("prob_Met_hi_EC", "mean"),
        ).reset_index()
        summary["analysis_set"] = changed_only
        for col in score_cols:
            vals = g[col].mean().reset_index(name=f"mean_{col}")
            summary = summary.merge(vals, on=["scenario", "metEClabel"], how="left")
        rows.append(summary)
    summary = pd.concat(rows, ignore_index=True)

    baseline = summary[summary["scenario"] == "baseline"].copy()
    delta = summary.merge(
        baseline.drop(columns=["scenario"]),
        on=["analysis_set", "metEClabel"],
        how="left",
        suffixes=("", "_baseline"),
    )
    for col in ["freq_pred_Met_hi_EC", "mean_prob_Met_hi_EC"] + [f"mean_{c}" for c in score_cols]:
        delta[f"delta_{col}"] = delta[col] - delta[f"{col}_baseline"]
    return delta


def summarise_samples(cell_df):
    score_cols = list(SCORE_GENESETS.keys())
    g = cell_df.groupby(["scenario", "sample_id", "metEClabel"], observed=True)
    sample = g.agg(
        n_ec=("cell_id", "size"),
        freq_pred_Met_hi_EC=("predicted_ec_state", lambda s: float((s == "Met_hi_EC").mean())),
        mean_prob_Met_hi_EC=("prob_Met_hi_EC", "mean"),
    ).reset_index()
    for col in score_cols:
        vals = g[col].mean().reset_index(name=col)
        sample = sample.merge(vals, on=["scenario", "sample_id", "metEClabel"], how="left")
    return sample


def summarise_ec_state_frequencies(cell_df):
    rows = []
    for analysis_set, df in [("all_ecs", cell_df), ("mcc_changed_ecs", cell_df[cell_df["mcc_changed_vs_baseline"]])]:
        if df.empty:
            continue
        total = df.groupby("scenario", observed=True)["cell_id"].size().rename("n_ec_total").reset_index()
        counts = (
            df.groupby(["scenario", "predicted_ec_state"], observed=True)["cell_id"]
            .size()
            .rename("n_ec_state")
            .reset_index()
        )
        out = counts.merge(total, on="scenario", how="left")
        out["freq_ec_state"] = out["n_ec_state"] / out["n_ec_total"]
        out["analysis_set"] = analysis_set
        rows.append(out)
    return pd.concat(rows, ignore_index=True)


def plot_ec_state_frequency_pairwise(ec_freq_df, outdir):
    sns.set_theme(style="whitegrid", context="talk")
    pdfs = []
    df = ec_freq_df[
        (ec_freq_df["analysis_set"] == "all_ecs")
        & ec_freq_df["scenario"].isin(
            ["baseline", "remove_Met_hi_Fib", "replace_Met_hi_Fib_with_Other_Fib"]
        )
    ].copy()
    baseline = df[df["scenario"] == "baseline"].copy()
    baseline_remove = pd.concat([baseline, df[df["scenario"] == "remove_Met_hi_Fib"]], ignore_index=True)
    baseline_remove["comparison"] = "Baseline vs remove Met_hi_Fib"
    baseline_replace = pd.concat(
        [baseline, df[df["scenario"] == "replace_Met_hi_Fib_with_Other_Fib"]],
        ignore_index=True,
    )
    baseline_replace["comparison"] = "Baseline vs replace Met_hi_Fib"
    pair_df = pd.concat([baseline_remove, baseline_replace], ignore_index=True)
    pair_df["scenario_label"] = pair_df["scenario"].map(
        {
            "baseline": "baseline",
            "remove_Met_hi_Fib": "remove",
            "replace_Met_hi_Fib_with_Other_Fib": "replace",
        }
    )
    pair_df["predicted_ec_state"] = pd.Categorical(pair_df["predicted_ec_state"], EC_LABELS)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), sharey=True)
    for ax, comparison in zip(axes, ["Baseline vs remove Met_hi_Fib", "Baseline vs replace Met_hi_Fib"]):
        cur = pair_df[pair_df["comparison"] == comparison]
        sns.barplot(
            data=cur,
            x="predicted_ec_state",
            y="freq_ec_state",
            hue="scenario_label",
            ax=ax,
        )
        for container in ax.containers:
            ax.bar_label(container, fmt="%.2f", fontsize=9, padding=3)
        ax.set_title(comparison)
        ax.set_xlabel("Predicted EC state among all ECs")
        ax.set_ylabel("Frequency among ECs")
        ax.set_ylim(0, 1.08)
    fig.tight_layout()
    p = outdir / "ec_state_frequency_among_all_ecs_pairwise.pdf"
    fig.savefig(p)
    fig.savefig(outdir / "ec_state_frequency_among_all_ecs_pairwise.png", dpi=220)
    pdfs.append(p)
    plt.close(fig)
    return pdfs


def plot_outputs(cell_df, summary_df, outdir):
    sns.set_theme(style="whitegrid", context="talk")
    pdfs = []

    score_long = cell_df.melt(
        id_vars=["scenario", "metEClabel", "cell_id"],
        value_vars=list(SCORE_GENESETS.keys()),
        var_name="score",
        value_name="value",
    )
    fig, ax = plt.subplots(figsize=(12, 6))
    sns.barplot(
        data=score_long,
        x="score",
        y="value",
        hue="scenario",
        errorbar=("ci", 95),
        ax=ax,
    )
    ax.set_xlabel("Score")
    ax.set_ylabel("Mean log1p generated Xmic score")
    ax.set_title("Generated metabolic scores across ECs")
    ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    p = outdir / "ec_metabolic_scores_all_ecs.pdf"
    fig.savefig(p)
    fig.savefig(outdir / "ec_metabolic_scores_all_ecs.png", dpi=180)
    pdfs.append(p)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(12, 6))
    sns.barplot(
        data=score_long,
        x="score",
        y="value",
        hue="metEClabel",
        errorbar=("ci", 95),
        ax=ax,
    )
    ax.set_xlabel("Score")
    ax.set_ylabel("Mean log1p generated Xmic score")
    ax.set_title("Generated metabolic scores by original EC population")
    ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    p = outdir / "ec_metabolic_scores_by_original_ec_population.pdf"
    fig.savefig(p)
    fig.savefig(outdir / "ec_metabolic_scores_by_original_ec_population.png", dpi=180)
    pdfs.append(p)
    plt.close(fig)

    return pdfs


def plot_tissue_maps(cell_df, tissues, outdir):
    sns.set_theme(style="white", context="talk")
    pdfs = []
    scenario_order = ["baseline", "remove_Met_hi_Fib", "replace_Met_hi_Fib_with_Other_Fib"]
    scenario_titles = {
        "baseline": "Baseline",
        "remove_Met_hi_Fib": "Remove Met_hi_Fib",
        "replace_Met_hi_Fib_with_Other_Fib": "Replace Met_hi_Fib with Other_Fib",
    }
    palette = {
        "Met_hi_Fib": "#c03a2b",
        "Other_Fib": "#7a7a7a",
        "Met_hi_EC": "#1f77b4",
        "Met_int_EC": "#2ca02c",
    }

    fig, axes = plt.subplots(1, 3, figsize=(26, 8), sharex=True, sharey=True)
    for ax, scenario in zip(axes, scenario_order):
        tissue = tissues[scenario]
        plot_df = tissue.obs[["combined_met_label", "x_centroid", "y_centroid"]].copy()
        plot_df["combined_met_label"] = plot_df["combined_met_label"].astype(str)
        for label, group in plot_df.groupby("combined_met_label", observed=True):
            ax.scatter(
                group["x_centroid"],
                group["y_centroid"],
                s=3,
                c=palette.get(label, "#333333"),
                label=label,
                linewidths=0,
                alpha=0.75,
                rasterized=True,
            )
        ax.set_title(scenario_titles[scenario])
        ax.set_aspect("equal", adjustable="box")
        ax.invert_yaxis()
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="center right", frameon=False)
    fig.tight_layout(rect=[0, 0, 0.9, 1])
    p = outdir / "tissue_maps_celltype_perturbations.pdf"
    fig.savefig(p)
    fig.savefig(outdir / "tissue_maps_celltype_perturbations.png", dpi=260)
    pdfs.append(p)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(26, 8), sharex=True, sharey=True)
    for ax, scenario in zip(axes, scenario_order):
        tissue = tissues[scenario]
        ec = cell_df[cell_df["scenario"] == scenario].merge(
            tissue.obs[["cell_id", "x_centroid", "y_centroid"]].reset_index(drop=True),
            on="cell_id",
            how="left",
        )
        base = tissue.obs[["x_centroid", "y_centroid"]]
        ax.scatter(
            base["x_centroid"],
            base["y_centroid"],
            s=1,
            c="#dddddd",
            linewidths=0,
            alpha=0.28,
            rasterized=True,
        )
        sca = ax.scatter(
            ec["x_centroid"],
            ec["y_centroid"],
            s=7,
            c=ec["prob_Met_hi_EC"],
            cmap="viridis",
            vmin=0,
            vmax=1,
            linewidths=0,
            alpha=0.92,
            rasterized=True,
        )
        ax.set_title(scenario_titles[scenario])
        ax.set_aspect("equal", adjustable="box")
        ax.invert_yaxis()
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    cbar = fig.colorbar(sca, ax=axes, fraction=0.025, pad=0.02)
    cbar.set_label("Generated EC probability: Met_hi_EC")
    fig.tight_layout()
    p = outdir / "tissue_maps_ec_met_hi_probability.pdf"
    fig.savefig(p)
    fig.savefig(outdir / "tissue_maps_ec_met_hi_probability.png", dpi=260)
    pdfs.append(p)
    plt.close(fig)

    for label in ["Met_hi_EC", "Met_int_EC"]:
        fig, axes = plt.subplots(1, 3, figsize=(26, 8), sharex=True, sharey=True)
        for ax, scenario in zip(axes, scenario_order):
            tissue = tissues[scenario]
            ec = cell_df[(cell_df["scenario"] == scenario) & (cell_df["metEClabel"] == label)].merge(
                tissue.obs[["cell_id", "x_centroid", "y_centroid"]].reset_index(drop=True),
                on="cell_id",
                how="left",
            )
            ax.scatter(
                tissue.obs["x_centroid"],
                tissue.obs["y_centroid"],
                s=1,
                c="#e5e5e5",
                linewidths=0,
                alpha=0.22,
                rasterized=True,
            )
            sca = ax.scatter(
                ec["x_centroid"],
                ec["y_centroid"],
                s=10,
                c=ec["prob_Met_hi_EC"],
                cmap="viridis",
                vmin=0,
                vmax=1,
                linewidths=0,
                alpha=0.95,
                rasterized=True,
            )
            ax.set_title(f"{scenario_titles[scenario]}: original {label}")
            ax.set_aspect("equal", adjustable="box")
            ax.invert_yaxis()
            ax.set_xlabel("x")
            ax.set_ylabel("y")
        cbar = fig.colorbar(sca, ax=axes, fraction=0.025, pad=0.02)
        cbar.set_label("Generated EC probability: Met_hi_EC")
        fig.tight_layout()
        safe_label = label.lower()
        p = outdir / f"tissue_maps_{safe_label}_met_hi_probability.pdf"
        fig.savefig(p)
        fig.savefig(outdir / f"tissue_maps_{safe_label}_met_hi_probability.png", dpi=260)
        pdfs.append(p)
        plt.close(fig)

    return pdfs


def write_html(outdir, summary_df, sample_df, ec_freq_df, used_genes, classifier_report, pdfs, command_args):
    html_path = outdir / "mintflow_met_hi_fib_perturbation_ec_report.html"
    summary_html = summary_df.to_html(index=False, float_format=lambda x: f"{x:.4g}")
    sample_html = sample_df.to_html(index=False, float_format=lambda x: f"{x:.4g}")
    ec_freq_html = ec_freq_df.to_html(index=False, float_format=lambda x: f"{x:.4g}")
    report_html = pd.DataFrame(classifier_report).T.to_html(float_format=lambda x: f"{x:.4g}")
    genes_html = "<ul>" + "".join(
        f"<li><b>{html.escape(k)}</b>: {html.escape(', '.join(v) if v else 'no genes found in panel')}</li>"
        for k, v in used_genes.items()
    ) + "</ul>"
    pdf_html = "<ul>" + "".join(
        f"<li>{html.escape(p.name)}</li>" for p in pdfs
    ) + "</ul>"
    pngs = [
        "ec_state_frequency_among_all_ecs_pairwise.png",
        "ec_metabolic_scores_all_ecs.png",
        "ec_metabolic_scores_by_original_ec_population.png",
        "tissue_maps_celltype_perturbations.png",
        "tissue_maps_ec_met_hi_probability.png",
        "tissue_maps_met_hi_ec_met_hi_probability.png",
        "tissue_maps_met_int_ec_met_hi_probability.png",
    ]
    imgs = "".join(f"<h3>{html.escape(p)}</h3><img src='{html.escape(p)}' style='max-width: 100%;'>" for p in pngs)
    content = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>MintFlow Met_hi_Fib Perturbation EC Analysis</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 32px; color: #222; line-height: 1.45; }}
    h1, h2, h3 {{ color: #111; }}
    table {{ border-collapse: collapse; margin: 16px 0; font-size: 13px; }}
    th, td {{ border: 1px solid #ccc; padding: 5px 7px; text-align: right; }}
    th:first-child, td:first-child {{ text-align: left; }}
    code {{ background: #f4f4f4; padding: 1px 4px; }}
    .note {{ background: #f7f7f7; border-left: 4px solid #777; padding: 12px 16px; }}
  </style>
</head>
<body>
  <h1>MintFlow Met_hi_Fib in silico perturbation: EC response</h1>
  <div class="note">
    <p><b>Input model:</b> {html.escape(command_args.checkpoint)}</p>
    <p><b>Input AnnData:</b> {html.escape(command_args.h5ad)}</p>
    <p><b>Generated expression:</b> average of {command_args.num_realisations} MintFlow realisations using <code>MintFLow_Generated_Xmic</code>, following the MintFlow in silico perturbation tutorial.</p>
  </div>

  <h2>What Was Perturbed</h2>
  <p>Three tissues were compared: baseline, Met_hi_Fib removed, and Met_hi_Fib relabelled as Other_Fib. For each tissue a fresh spatial neighbour graph was computed before calling <code>mintflow.generate_insilico_ST_data</code>.</p>

  <h2>How EC Frequencies Were Estimated</h2>
    <p>MintFlow generation keeps supplied cell type labels fixed. To estimate whether generated EC expression resembles Met_hi_EC or Met_int_EC, a logistic-regression classifier was trained on the baseline MintFlow-generated Xmic expression of ECs using the original transferred labels <code>metEClabel</code>. The same classifier was then applied to the remove and replace perturbations. The primary frequency table reports predicted EC-state composition among all ECs, so baseline starts from the observed EC mixture rather than per-class recall.</p>
  <h3>Classifier Training Performance on Baseline Generated ECs</h3>
  {report_html}

  <h2>Metabolic Scores</h2>
  <p>Scores are mean log1p generated Xmic expression across marker genes present in the Xenium panel. Marker mappings were adapted from the IMC score columns requested in the prompt.</p>
  {genes_html}

  <h2>Main Summary</h2>
  <h3>EC-state frequencies among all ECs</h3>
  {ec_freq_html}
  <h3>Per-original-EC-population classifier shifts</h3>
  {summary_html}

  <h2>Sample-Level Summary</h2>
  {sample_html}

  <h2>Plots</h2>
  <p>PDF exports written separately:</p>
  {pdf_html}
  {imgs}

  <h2>Run Parameters</h2>
  <pre>{html.escape(json.dumps(vars(command_args), indent=2))}</pre>
</body>
</html>
"""
    html_path.write_text(content)
    return html_path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--num-realisations", type=int, default=5)
    parser.add_argument("--n-neighs", type=int, default=5)
    parser.add_argument("--reuse-generations", action="store_true")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    ensure_dir(outdir)

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    checkpoint = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    checkpoint["model"].to(device)
    checkpoint["model"].eval()

    adata = ad.read_h5ad(args.h5ad)
    adata.obs["cell_id"] = adata.obs["cell_id"].astype(str)
    adata.obs["combined_met_label"] = adata.obs["combined_met_label"].astype(str)
    adata.obs["metEClabel"] = adata.obs["metEClabel"].astype(str)

    tissues = {}
    tissues["baseline"] = adata.copy()
    tissues["remove_Met_hi_Fib"] = adata[adata.obs["combined_met_label"] != "Met_hi_Fib"].copy()
    tissues["replace_Met_hi_Fib_with_Other_Fib"] = adata.copy()
    tissues["replace_Met_hi_Fib_with_Other_Fib"].obs["combined_met_label"] = tissues[
        "replace_Met_hi_Fib_with_Other_Fib"
    ].obs["combined_met_label"].replace({"Met_hi_Fib": "Other_Fib"})

    results = {}
    for name, tissue in tissues.items():
        make_spatial_graph(tissue, n_neighs=args.n_neighs)
        generation_path = outdir / f"generation_{name}.pkl"
        if args.reuse_generations and generation_path.exists():
            print(f"Reusing MintFlow generation for {name}: {generation_path}")
            with open(generation_path, "rb") as handle:
                results[name] = pickle.load(handle)
        else:
            print(f"Generating MintFlow expression for {name}: {tissue.n_obs} cells")
            results[name] = mintflow.generate_insilico_ST_data(
                adata=tissue,
                obskey_celltype="combined_met_label",
                obspkey_neighbourhood_graph="spatial_connectivities",
                device=device,
                batch_index_trainingdata=0,
                num_generated_realisations=args.num_realisations,
                model=checkpoint["model"],
                data_mintflow=checkpoint["data_mintflow"],
                dict_all4_configs=checkpoint["dict_all4_configs"],
                estimate_spatial_sizefactors_on_sections=[0],
            )
            with open(generation_path, "wb") as handle:
                pickle.dump(results[name], handle)

    baseline_xmic = avg_generated_xmic(results["baseline"])
    clf, classifier_report = train_ec_classifier_from_generated(tissues["baseline"], baseline_xmic)

    cell_tables = []
    used_genes_final = None
    for name, tissue in tissues.items():
        Xmic = avg_generated_xmic(results[name])
        cell_df, used_genes = analyze_scenario(
            name=name,
            adata_scenario=tissue,
            Xmic=Xmic,
            clf=clf,
            adata_base=tissues["baseline"],
            result_base=results["baseline"],
            result_scenario=results[name],
        )
        cell_tables.append(cell_df)
        used_genes_final = used_genes

    cell_df = pd.concat(cell_tables, ignore_index=True)
    summary_df = summarise_cells(cell_df)
    sample_df = summarise_samples(cell_df)
    ec_freq_df = summarise_ec_state_frequencies(cell_df)

    cell_df.to_csv(outdir / "ec_cell_level_perturbation_results.csv", index=False)
    summary_df.to_csv(outdir / "ec_summary_perturbation_results.csv", index=False)
    sample_df.to_csv(outdir / "ec_sample_level_perturbation_results.csv", index=False)
    ec_freq_df.to_csv(outdir / "ec_state_frequency_among_all_ecs.csv", index=False)
    Path(outdir / "score_genes_used.json").write_text(json.dumps(used_genes_final, indent=2))
    Path(outdir / "ec_classifier_report.json").write_text(json.dumps(classifier_report, indent=2))

    pdfs = plot_outputs(cell_df, summary_df, outdir)
    pdfs.extend(plot_ec_state_frequency_pairwise(ec_freq_df, outdir))
    pdfs.extend(plot_tissue_maps(cell_df, tissues, outdir))
    html_path = write_html(outdir, summary_df, sample_df, ec_freq_df, used_genes_final, classifier_report, pdfs, args)
    print(f"Finished. HTML report: {html_path}")


if __name__ == "__main__":
    main()
