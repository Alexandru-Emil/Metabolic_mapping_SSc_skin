
import os
from pathlib import Path
def project_path(relative):
    return Path(os.environ.get("METABOLIC_INPUT_DIR", "inputs")) / relative

#!/usr/bin/env python

from pathlib import Path
import html
import pickle

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import sparse, stats
from scipy.spatial import cKDTree
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    confusion_matrix,
    roc_auc_score,
)
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


WORKDIR = Path(project_path('external/mintflow_xenium_training'))
H5AD = WORKDIR / "xenium_sfe_09032026_combined_fib_ec_metlabel.h5ad"
BASELINE_GENERATION = WORKDIR / "perturb_met_hi_fib_ec" / "generation_baseline.pkl"
OUT = WORKDIR / "cross_modal_ec_validation_plots"

EC_LABELS = ["Met_hi_EC", "Met_int_EC"]
LABEL_COLORS = {"Met_hi_EC": "#2b7bba", "Met_int_EC": "#55a868"}
SCORE_GENESETS = {
    "glycolysis": ["SLC2A1", "GLUT1", "HK1", "PKM", "PKM2", "LDHA"],
    "tca_oxphos": ["OGDH", "ATP5F1A", "ATP5A", "SDHA"],
    "hif1a": ["HIF1A"],
    "nox4": ["NOX4"],
    "vascular_activation": ["KDR", "FLT1", "VWF", "PECAM1", "ENG", "ICAM1", "SELE", "VCAM1", "ANGPT2"],
}
PAPER_RC = {
    "axes.titlesize": 12,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 8,
    "legend.title_fontsize": 9,
}


def set_theme(style="whitegrid"):
    sns.set_theme(style=style, context="paper", font_scale=1.15, rc=PAPER_RC)


def savefig(fig, stem, dpi=300):
    pdf = OUT / f"{stem}.pdf"
    png = OUT / f"{stem}.png"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return pdf, png


def as_dense(X):
    return X.toarray() if sparse.issparse(X) else np.asarray(X)


def lognorm_counts(X, scale_factor=1e4):
    X = X.copy()
    if sparse.issparse(X):
        lib = np.asarray(X.sum(axis=1)).ravel()
        lib[lib == 0] = 1
        X = sparse.diags(scale_factor / lib).dot(X)
        X = X.tocsr()
        X.data = np.log1p(X.data)
        return X
    X = np.asarray(X, dtype=float)
    lib = X.sum(axis=1)
    lib[lib == 0] = 1
    return np.log1p(X / lib[:, None] * scale_factor)


def as_log1p_nonnegative(X):
    if sparse.issparse(X):
        X = X.copy()
        X.data[X.data < 0] = 0
        X.data = np.log1p(X.data)
        return X
    return np.log1p(np.clip(np.asarray(X), 0, None))


def avg_generated_xmic(result):
    return np.stack(
        [r["MintFLow_Generated_Xmic"] for r in result["list_generated_realisations_ie_expressions"]],
        axis=0,
    ).mean(axis=0)


def bh_adjust(pvals):
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    ok = np.isfinite(p)
    p_ok = p[ok]
    n = len(p_ok)
    if n == 0:
        return out
    order = np.argsort(p_ok)
    ranks = np.arange(1, n + 1)
    adj_sorted = p_ok[order] * n / ranks
    adj_sorted = np.minimum.accumulate(adj_sorted[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.minimum(adj_sorted, 1.0)
    out[ok] = adj
    out[~ok] = 1.0
    return out


def gene_lookup(var):
    lookup = {}
    for idx, row in var.iterrows():
        for val in [idx, row.get("gene_id", ""), row.get("gene_symbol", "")]:
            if pd.notna(val):
                lookup[str(val).upper()] = idx
    return lookup


def genes_present(var, genes):
    lookup = gene_lookup(var)
    seen = []
    for gene in genes:
        hit = lookup.get(gene.upper())
        if hit is not None and hit not in seen:
            seen.append(hit)
    return seen


def load_data():
    adata = ad.read_h5ad(H5AD)
    adata.obs["cell_id"] = adata.obs["cell_id"].astype(str)
    adata.obs["combined_met_label"] = adata.obs["combined_met_label"].astype(str)
    adata.obs["metEClabel"] = adata.obs["metEClabel"].astype(str)
    adata.obs["sample_id"] = adata.obs["sample_id"].astype(str)
    ec = adata[adata.obs["metEClabel"].isin(EC_LABELS)].copy()
    ec.obs["metEClabel"] = pd.Categorical(ec.obs["metEClabel"], EC_LABELS, ordered=True)

    sample_ids = sorted(ec.obs["sample_id"].unique())
    sample_map = {sid: f"S{i + 1}" for i, sid in enumerate(sample_ids)}
    ec.obs["sample_short"] = ec.obs["sample_id"].map(sample_map)
    adata.obs["sample_short"] = adata.obs["sample_id"].map(sample_map)
    pd.DataFrame({"sample_id": list(sample_map), "sample_short": list(sample_map.values())}).to_csv(
        OUT / "sample_id_short_label_map.csv", index=False
    )

    X_log = lognorm_counts(ec.X)
    X_log_dense = as_dense(X_log)
    return adata, ec, X_log, X_log_dense


def panel_a_schematic():
    set_theme("white")
    fig, ax = plt.subplots(figsize=(10, 4.2))
    ax.axis("off")
    boxes = [
        (0.04, 0.55, 0.22, 0.25, "IMC protein data\nmetabolic markers"),
        (0.39, 0.55, 0.23, 0.25, "IMC-defined EC labels\nMet_hi_EC / Met_int_EC"),
        (0.74, 0.55, 0.22, 0.25, "Matched Xenium ECs\nobserved RNA profiles"),
        (0.39, 0.15, 0.23, 0.22, "Shared cell KEY\ncell-level matching"),
        (0.74, 0.15, 0.22, 0.22, "RNA-space validation\nDE, scores, classifier"),
    ]
    for x, y, w, h, label in boxes:
        rect = plt.Rectangle((x, y), w, h, facecolor="#f6f6f6", edgecolor="#333333", linewidth=1.4)
        ax.add_patch(rect)
        ax.text(x + w / 2, y + h / 2, label, ha="center", va="center", fontsize=11)
    arrows = [
        ((0.26, 0.68), (0.39, 0.68), "define"),
        ((0.62, 0.68), (0.74, 0.68), "transfer"),
        ((0.505, 0.55), (0.505, 0.37), "matched by KEY"),
        ((0.62, 0.26), (0.74, 0.26), "validate"),
        ((0.85, 0.55), (0.85, 0.37), ""),
    ]
    for start, end, label in arrows:
        ax.annotate("", xy=end, xytext=start, arrowprops=dict(arrowstyle="->", linewidth=1.6, color="#333333"))
        if label:
            ax.text((start[0] + end[0]) / 2, (start[1] + end[1]) / 2 + 0.04, label, ha="center", fontsize=9)
    ax.set_title("IMC-to-Xenium EC metabolic label transfer and RNA-space validation")
    return savefig(fig, "panel_A_imc_to_xenium_label_transfer_schematic")


def panel_b_embedding(ec, X_log_dense):
    set_theme("whitegrid")
    n_pcs = min(30, X_log_dense.shape[0] - 1, X_log_dense.shape[1])
    pca = PCA(n_components=n_pcs, random_state=0)
    pcs = pca.fit_transform(X_log_dense)
    ec.obs["PC1"] = pcs[:, 0]
    ec.obs["PC2"] = pcs[:, 1]
    ec.obsm["X_pca"] = pcs
    try:
        sc.pp.neighbors(ec, use_rep="X_pca", n_neighbors=18, random_state=0)
        sc.tl.umap(ec, random_state=0, min_dist=0.35)
        has_umap = True
    except Exception as exc:
        print(f"UMAP failed; PCA will still be plotted. Reason: {exc}")
        has_umap = False

    fig, axes = plt.subplots(1, 2 if has_umap else 1, figsize=(10.5 if has_umap else 5.4, 4.6))
    axes = np.atleast_1d(axes)
    sns.scatterplot(
        data=ec.obs,
        x="PC1",
        y="PC2",
        hue="metEClabel",
        palette=LABEL_COLORS,
        s=34,
        linewidth=0,
        alpha=0.85,
        ax=axes[0],
    )
    axes[0].set_title("Observed EC PCA")
    axes[0].set_xlabel(f"PC1 ({100 * pca.explained_variance_ratio_[0]:.1f}%)")
    axes[0].set_ylabel(f"PC2 ({100 * pca.explained_variance_ratio_[1]:.1f}%)")
    axes[0].legend(title="IMC EC label", frameon=False)
    if has_umap:
        um = ec.obsm["X_umap"]
        plot_df = ec.obs.copy()
        plot_df["UMAP1"] = um[:, 0]
        plot_df["UMAP2"] = um[:, 1]
        sns.scatterplot(
            data=plot_df,
            x="UMAP1",
            y="UMAP2",
            hue="metEClabel",
            palette=LABEL_COLORS,
            s=34,
            linewidth=0,
            alpha=0.85,
            ax=axes[1],
        )
        axes[1].set_title("Observed EC UMAP")
        axes[1].legend(title="IMC EC label", frameon=False)
    fig.tight_layout()
    return savefig(fig, "panel_B_observed_xenium_ec_embedding_by_imc_label")


def pseudobulk_profiles(ec):
    X_counts = as_dense(ec.X)
    genes = ec.var["gene_symbol"].astype(str).to_numpy() if "gene_symbol" in ec.var.columns else ec.var_names.astype(str)
    rows = []
    mat = []
    for (sample, state), idx in ec.obs.groupby(["sample_id", "metEClabel"], observed=True).indices.items():
        counts = X_counts[list(idx), :].sum(axis=0)
        lib = counts.sum() or 1
        logcpm = np.log1p(counts / lib * 1e6)
        rows.append({"sample_id": sample, "sample_short": ec.obs.loc[ec.obs["sample_id"] == sample, "sample_short"].iloc[0], "metEClabel": state, "n_cells": len(idx)})
        mat.append(logcpm)
    meta = pd.DataFrame(rows)
    pb = pd.DataFrame(np.vstack(mat), columns=genes)
    return meta, pb


def panel_c_pseudobulk_pca(meta, pb):
    set_theme("whitegrid")
    X = StandardScaler().fit_transform(pb.to_numpy())
    pcs = PCA(n_components=2, random_state=0).fit_transform(X)
    pca = PCA(n_components=2, random_state=0).fit(X)
    plot_df = meta.copy()
    plot_df["PC1"] = pcs[:, 0]
    plot_df["PC2"] = pcs[:, 1]
    plot_df.to_csv(OUT / "panel_C_pseudobulk_sample_state_pca_coordinates.csv", index=False)

    fig, ax = plt.subplots(figsize=(6.6, 5.2))
    sns.scatterplot(
        data=plot_df,
        x="PC1",
        y="PC2",
        hue="metEClabel",
        style="sample_short",
        palette=LABEL_COLORS,
        s=95,
        ax=ax,
    )
    for _, row in plot_df.iterrows():
        ax.text(row["PC1"] + 0.04, row["PC2"] + 0.04, row["sample_short"], fontsize=8)
    ax.set_xlabel(f"PC1 ({100 * pca.explained_variance_ratio_[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({100 * pca.explained_variance_ratio_[1]:.1f}%)")
    ax.set_title("Pseudobulk PCA: sample x IMC EC state")
    ax.legend(title="", frameon=False, loc="center left", bbox_to_anchor=(1.01, 0.5))
    fig.tight_layout()
    return savefig(fig, "panel_C_pseudobulk_pca_sample_by_ec_state")


def differential_expression(ec, meta, pb):
    genes = pb.columns.to_numpy()
    X_counts = as_dense(ec.X)
    rows = []
    for gene_idx, gene in enumerate(genes):
        hi_vals = []
        int_vals = []
        for sample in sorted(meta["sample_id"].unique()):
            hi = meta.index[(meta["sample_id"] == sample) & (meta["metEClabel"] == "Met_hi_EC")]
            mid = meta.index[(meta["sample_id"] == sample) & (meta["metEClabel"] == "Met_int_EC")]
            if len(hi) and len(mid):
                hi_vals.append(float(pb.iloc[hi[0], gene_idx]))
                int_vals.append(float(pb.iloc[mid[0], gene_idx]))
        hi_vals = np.array(hi_vals)
        int_vals = np.array(int_vals)
        lfc = float(np.mean(hi_vals - int_vals)) if len(hi_vals) else np.nan
        if len(hi_vals) >= 3 and np.nanstd(hi_vals - int_vals) > 0:
            pval = float(stats.ttest_rel(hi_vals, int_vals, nan_policy="omit").pvalue)
        else:
            pval = 1.0
        hi_cells = ec.obs["metEClabel"].astype(str).to_numpy() == "Met_hi_EC"
        int_cells = ec.obs["metEClabel"].astype(str).to_numpy() == "Met_int_EC"
        counts = X_counts[:, gene_idx]
        rows.append(
            {
                "gene": gene,
                "logFC_Met_hi_vs_Met_int": lfc,
                "pvalue": pval,
                "mean_logCPM_Met_hi_EC": float(np.mean(hi_vals)) if len(hi_vals) else np.nan,
                "mean_logCPM_Met_int_EC": float(np.mean(int_vals)) if len(int_vals) else np.nan,
                "frac_expr_Met_hi_EC": float(np.mean(counts[hi_cells] > 0)),
                "frac_expr_Met_int_EC": float(np.mean(counts[int_cells] > 0)),
                "n_paired_samples": len(hi_vals),
            }
        )
    de = pd.DataFrame(rows)
    de["padj"] = bh_adjust(de["pvalue"])
    de["neglog10_padj"] = -np.log10(np.maximum(de["padj"].to_numpy(), 1e-300))
    de = de.sort_values(["padj", "pvalue", "logFC_Met_hi_vs_Met_int"], ascending=[True, True, False])
    de.to_csv(OUT / "panel_D_observed_ec_pseudobulk_DE_Met_hi_vs_Met_int.csv", index=False)
    de.nlargest(30, "logFC_Met_hi_vs_Met_int").to_csv(OUT / "panel_D_top_Met_hi_EC_enriched_genes.csv", index=False)
    de.nsmallest(30, "logFC_Met_hi_vs_Met_int").to_csv(OUT / "panel_D_top_Met_int_EC_enriched_genes.csv", index=False)
    return de


def panel_d_volcano(de):
    set_theme("whitegrid")
    plot_df = de.copy()
    plot_df["direction"] = "not significant"
    plot_df.loc[(plot_df["padj"] < 0.1) & (plot_df["logFC_Met_hi_vs_Met_int"] > 0), "direction"] = "Met_hi_EC enriched"
    plot_df.loc[(plot_df["padj"] < 0.1) & (plot_df["logFC_Met_hi_vs_Met_int"] < 0), "direction"] = "Met_int_EC enriched"
    colors = {"not significant": "#bdbdbd", "Met_hi_EC enriched": "#2b7bba", "Met_int_EC enriched": "#55a868"}

    fig, ax = plt.subplots(figsize=(6.7, 5.4))
    sns.scatterplot(
        data=plot_df,
        x="logFC_Met_hi_vs_Met_int",
        y="neglog10_padj",
        hue="direction",
        palette=colors,
        s=14,
        linewidth=0,
        alpha=0.75,
        ax=ax,
    )
    top = pd.concat([
        de.nlargest(6, "logFC_Met_hi_vs_Met_int"),
        de.nsmallest(6, "logFC_Met_hi_vs_Met_int"),
    ])
    for _, row in top.iterrows():
        ax.text(row["logFC_Met_hi_vs_Met_int"], row["neglog10_padj"] + 0.05, row["gene"], fontsize=7, ha="center")
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Pseudobulk logFC: Met_hi_EC - Met_int_EC")
    ax.set_ylabel("-log10 adjusted P value")
    ax.set_title("Observed Xenium EC pseudobulk DE")
    ax.legend(title="", frameon=False, loc="center left", bbox_to_anchor=(1.01, 0.5))
    fig.tight_layout()
    return savefig(fig, "panel_D_volcano_observed_ec_met_hi_vs_met_int")


def top_de_genes(de, n=8):
    hi = de.sort_values("logFC_Met_hi_vs_Met_int", ascending=False).head(n)["gene"].tolist()
    mid = de.sort_values("logFC_Met_hi_vs_Met_int", ascending=True).head(n)["gene"].tolist()
    return hi + [g for g in mid if g not in hi]


def panel_e_top_gene_heatmap_dotplot(ec, X_log_dense, de):
    genes = top_de_genes(de, n=8)
    gene_to_idx = {str(g): i for i, g in enumerate(ec.var["gene_symbol"].astype(str))}
    idx = [gene_to_idx[g] for g in genes if g in gene_to_idx]
    genes = [genes[i] for i, g in enumerate(genes) if g in gene_to_idx]

    rows = []
    for sample in sorted(ec.obs["sample_id"].unique()):
        for state in EC_LABELS:
            mask = (ec.obs["sample_id"].astype(str).to_numpy() == sample) & (ec.obs["metEClabel"].astype(str).to_numpy() == state)
            if mask.sum() == 0:
                continue
            vals = X_log_dense[mask][:, idx].mean(axis=0)
            row = {"sample_short": ec.obs.loc[ec.obs["sample_id"] == sample, "sample_short"].iloc[0], "metEClabel": state}
            row.update({g: vals[i] for i, g in enumerate(genes)})
            rows.append(row)
    heat = pd.DataFrame(rows)
    mat = heat[genes].to_numpy()
    z = (mat - mat.mean(axis=0)) / (mat.std(axis=0) + 1e-9)
    heat_z = pd.DataFrame(z, columns=genes)
    heat_z.index = heat["sample_short"] + "_" + heat["metEClabel"].astype(str)
    heat.to_csv(OUT / "panel_E_top_gene_pseudobulk_expression_by_sample_state.csv", index=False)

    dot_rows = []
    counts = as_dense(ec.X)
    for state in EC_LABELS:
        mask = ec.obs["metEClabel"].astype(str).to_numpy() == state
        for gene, j in zip(genes, idx):
            dot_rows.append(
                {
                    "gene": gene,
                    "metEClabel": state,
                    "mean_log_expression": float(X_log_dense[mask, j].mean()),
                    "fraction_expressing": float((counts[mask, j] > 0).mean()),
                }
            )
    dot = pd.DataFrame(dot_rows)
    dot.to_csv(OUT / "panel_E_top_gene_dotplot_values.csv", index=False)

    set_theme("white")
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 6.0), gridspec_kw={"width_ratios": [1.4, 1.0]})
    sns.heatmap(heat_z, cmap="vlag", center=0, linewidths=0.2, linecolor="white", ax=axes[0], cbar_kws={"label": "z-scored pseudobulk logCPM"})
    axes[0].set_title("Top DE genes by sample and EC state")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("")

    sca = axes[1].scatter(
        x=pd.Categorical(dot["metEClabel"], EC_LABELS).codes,
        y=pd.Categorical(dot["gene"], genes[::-1]).codes,
        s=40 + 260 * dot["fraction_expressing"],
        c=dot["mean_log_expression"],
        cmap="viridis",
        edgecolor="black",
        linewidth=0.3,
    )
    axes[1].set_xticks(range(len(EC_LABELS)))
    axes[1].set_xticklabels(EC_LABELS, rotation=20)
    axes[1].set_yticks(range(len(genes)))
    axes[1].set_yticklabels(genes[::-1])
    axes[1].set_title("Mean expression and fraction expressing")
    axes[1].set_xlabel("")
    axes[1].set_ylabel("")
    cbar = fig.colorbar(sca, ax=axes[1], fraction=0.05, pad=0.04)
    cbar.set_label("Mean log-normalized expression")
    fig.tight_layout()
    return savefig(fig, "panel_E_top_de_gene_heatmap_and_dotplot")


def compute_scores(ec, X_log_dense):
    score_df = ec.obs[["cell_id", "sample_id", "sample_short", "metEClabel", "x_centroid", "y_centroid"]].copy()
    used = {}
    for score, genes in SCORE_GENESETS.items():
        present = genes_present(ec.var, genes)
        used[score] = present
        if not present:
            score_df[score] = np.nan
            continue
        idx = [list(ec.var_names).index(g) for g in present]
        score_df[score] = X_log_dense[:, idx].mean(axis=1)
    pd.DataFrame(
        [{"score": k, "genes_present": ", ".join(v) if v else ""} for k, v in used.items()]
    ).to_csv(OUT / "panel_F_observed_rna_score_genes_used.csv", index=False)
    score_df.to_csv(OUT / "panel_F_observed_ec_rna_scores.csv", index=False)
    return score_df


def panel_f_scores(score_df):
    score_cols = list(SCORE_GENESETS)
    long = score_df.melt(
        id_vars=["cell_id", "sample_id", "sample_short", "metEClabel"],
        value_vars=score_cols,
        var_name="score",
        value_name="value",
    ).dropna()
    set_theme("whitegrid")
    fig, axes = plt.subplots(1, len(score_cols), figsize=(16, 4.4), sharey=False)
    for ax, score in zip(axes, score_cols):
        cur = long[long["score"] == score]
        sns.boxplot(
            data=cur,
            x="metEClabel",
            y="value",
            hue="metEClabel",
            palette=LABEL_COLORS,
            fliersize=0,
            width=0.55,
            ax=ax,
            legend=False,
        )
        sns.stripplot(data=cur, x="metEClabel", y="value", color="black", size=2, alpha=0.35, jitter=0.18, ax=ax)
        ax.set_title(score)
        ax.set_xlabel("")
        ax.set_ylabel("Observed RNA score")
        ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    return savefig(fig, "panel_F_observed_rna_metabolic_stress_scores_by_imc_ec_state")


def panel_g_sample_paired_scores(score_df):
    score_cols = list(SCORE_GENESETS)
    sample = (
        score_df.groupby(["sample_id", "sample_short", "metEClabel"], observed=True)[score_cols]
        .mean()
        .reset_index()
    )
    sample.to_csv(OUT / "panel_G_sample_level_observed_ec_rna_scores.csv", index=False)
    long = sample.melt(
        id_vars=["sample_id", "sample_short", "metEClabel"],
        value_vars=score_cols,
        var_name="score",
        value_name="value",
    )
    set_theme("whitegrid")
    fig, axes = plt.subplots(1, len(score_cols), figsize=(16, 4.4), sharey=False)
    for ax, score in zip(axes, score_cols):
        cur = long[long["score"] == score].copy()
        for _, g in cur.groupby("sample_id", observed=True):
            g = g.sort_values("metEClabel")
            ax.plot(g["metEClabel"], g["value"], color="#777777", alpha=0.65, linewidth=1)
        sns.stripplot(
            data=cur,
            x="metEClabel",
            y="value",
            hue="metEClabel",
            palette=LABEL_COLORS,
            size=6,
            jitter=0.05,
            edgecolor="white",
            linewidth=0.4,
            ax=ax,
            legend=False,
        )
        ax.set_title(score)
        ax.set_xlabel("")
        ax.set_ylabel("Sample mean RNA score")
        ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    return savefig(fig, "panel_G_per_sample_paired_observed_rna_score_differences")


def train_classifier(X, y):
    return make_pipeline(
        StandardScaler(with_mean=False),
        LogisticRegression(max_iter=5000, class_weight="balanced", solver="liblinear"),
    ).fit(X, y)


def panel_h_loso_classifier(ec, X_log):
    y = ec.obs["metEClabel"].astype(str).to_numpy()
    y_bin = (y == "Met_hi_EC").astype(int)
    samples = ec.obs["sample_id"].astype(str).to_numpy()
    rows = []
    all_true = []
    all_pred = []
    all_prob = []
    for sample in sorted(np.unique(samples)):
        test = samples == sample
        train = ~test
        if len(np.unique(y[train])) < 2 or len(np.unique(y[test])) < 2:
            continue
        clf = train_classifier(X_log[train], y[train])
        pred = clf.predict(X_log[test])
        classes = list(clf.named_steps["logisticregression"].classes_)
        prob = clf.predict_proba(X_log[test])[:, classes.index("Met_hi_EC")]
        rows.append(
            {
                "sample_id": sample,
                "sample_short": ec.obs.loc[ec.obs["sample_id"] == sample, "sample_short"].iloc[0],
                "n_test": int(test.sum()),
                "balanced_accuracy": balanced_accuracy_score(y[test], pred),
                "auroc": roc_auc_score(y_bin[test], prob),
                "average_precision": average_precision_score(y_bin[test], prob),
            }
        )
        all_true.extend(y[test])
        all_pred.extend(pred)
        all_prob.extend(prob)
    perf = pd.DataFrame(rows)
    perf.to_csv(OUT / "panel_H_leave_one_sample_out_classifier_metrics.csv", index=False)
    cm = confusion_matrix(all_true, all_pred, labels=EC_LABELS)
    pd.DataFrame(cm, index=EC_LABELS, columns=EC_LABELS).to_csv(OUT / "panel_H_leave_one_sample_out_confusion_matrix.csv")

    set_theme("whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.4))
    perf_long = perf.melt(
        id_vars=["sample_id", "sample_short"],
        value_vars=["balanced_accuracy", "auroc", "average_precision"],
        var_name="metric",
        value_name="value",
    )
    sns.stripplot(data=perf_long, x="metric", y="value", color="#333333", size=7, jitter=0.18, ax=axes[0])
    sns.pointplot(data=perf_long, x="metric", y="value", errorbar=("pi", 50), color="#2b7bba", join=False, ax=axes[0])
    axes[0].axhline(0.5, color="black", linestyle="--", linewidth=1)
    axes[0].set_ylim(0, 1.05)
    axes[0].set_xlabel("")
    axes[0].set_ylabel("Leave-one-sample-out performance")
    axes[0].tick_params(axis="x", rotation=25)
    axes[0].set_title("Observed RNA predicts IMC EC state")
    sns.heatmap(cm, annot=True, fmt="d", cmap="Blues", xticklabels=EC_LABELS, yticklabels=EC_LABELS, ax=axes[1], cbar=False)
    axes[1].set_xlabel("Predicted")
    axes[1].set_ylabel("Observed IMC label")
    axes[1].set_title("Aggregated confusion matrix")
    fig.tight_layout()
    return savefig(fig, "panel_H_leave_one_sample_out_observed_rna_classifier")


def panel_i_observed_vs_generated_classifier(ec, X_log, adata):
    y = ec.obs["metEClabel"].astype(str).to_numpy()
    obs_clf = train_classifier(X_log, y)
    obs_classes = list(obs_clf.named_steps["logisticregression"].classes_)
    prob_obs = obs_clf.predict_proba(X_log)[:, obs_classes.index("Met_hi_EC")]

    with BASELINE_GENERATION.open("rb") as handle:
        baseline_result = pickle.load(handle)
    Xmic = avg_generated_xmic(baseline_result)
    ec_mask_all = adata.obs["metEClabel"].isin(EC_LABELS).to_numpy()
    gen_clf = train_classifier(as_log1p_nonnegative(Xmic[ec_mask_all, :]), y)
    gen_classes = list(gen_clf.named_steps["logisticregression"].classes_)
    prob_gen = gen_clf.predict_proba(as_log1p_nonnegative(Xmic[ec_mask_all, :]))[:, gen_classes.index("Met_hi_EC")]

    concord = ec.obs[["cell_id", "metEClabel"]].copy()
    concord["observed_rna_classifier_prob_Met_hi_EC"] = prob_obs
    concord["baseline_generated_classifier_prob_Met_hi_EC"] = prob_gen
    concord.to_csv(OUT / "panel_I_observed_vs_generated_classifier_probabilities.csv", index=False)
    corr = np.corrcoef(prob_obs, prob_gen)[0, 1]

    set_theme("whitegrid")
    fig, ax = plt.subplots(figsize=(5.7, 5.2))
    sns.scatterplot(
        data=concord,
        x="observed_rna_classifier_prob_Met_hi_EC",
        y="baseline_generated_classifier_prob_Met_hi_EC",
        hue="metEClabel",
        palette=LABEL_COLORS,
        s=32,
        linewidth=0,
        alpha=0.75,
        ax=ax,
    )
    ax.plot([0, 1], [0, 1], linestyle="--", color="black", linewidth=1)
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlabel("Observed RNA classifier probability Met_hi_EC")
    ax.set_ylabel("Baseline generated classifier probability Met_hi_EC")
    ax.set_title(f"Observed vs generated classifier concordance\nPearson r = {corr:.2f}")
    ax.legend(title="IMC EC label", frameon=False, loc="lower right")
    fig.tight_layout()
    return savefig(fig, "panel_I_observed_vs_generated_ec_classifier_concordance")


def panel_j_spatial_scores(score_df, ec, X_log):
    y = ec.obs["metEClabel"].astype(str).to_numpy()
    clf = train_classifier(X_log, y)
    classes = list(clf.named_steps["logisticregression"].classes_)
    score_df = score_df.copy()
    score_df["observed_classifier_prob_Met_hi_EC"] = clf.predict_proba(X_log)[:, classes.index("Met_hi_EC")]
    score_df.to_csv(OUT / "panel_J_observed_ec_scores_and_classifier_probabilities.csv", index=False)

    set_theme("white")
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.4), sharex=True, sharey=True)
    for ax, col, title in [
        (axes[0], "glycolysis", "Observed EC glycolysis RNA score"),
        (axes[1], "observed_classifier_prob_Met_hi_EC", "Observed RNA classifier probability"),
    ]:
        sca = ax.scatter(
            score_df["x_centroid"],
            score_df["y_centroid"],
            c=score_df[col],
            cmap="viridis",
            s=18,
            linewidths=0.35,
            edgecolors=score_df["metEClabel"].astype(str).map({"Met_hi_EC": "#0b3d75", "Met_int_EC": "#1f5f2e"}),
            alpha=0.95,
        )
        cbar = fig.colorbar(sca, ax=ax, fraction=0.04, pad=0.02)
        cbar.set_label(col)
        ax.set_title(title)
        ax.set_aspect("equal", adjustable="box")
        ax.invert_yaxis()
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    fig.tight_layout()
    return savefig(fig, "panel_J_spatial_map_observed_ec_rna_scores")


def panel_k_proximity(adata, ec):
    fib = adata.obs[adata.obs["combined_met_label"] == "Met_hi_Fib"][["x_centroid", "y_centroid"]].to_numpy(float)
    ec_xy = ec.obs[["x_centroid", "y_centroid"]].to_numpy(float)
    dist = cKDTree(fib).query(ec_xy, k=1)[0]
    prox = ec.obs[["cell_id", "sample_id", "sample_short", "metEClabel"]].copy()
    prox["distance_to_nearest_Met_hi_Fib"] = dist
    bins = [0, 10, 20, 40, 80, 160, np.inf]
    labels = ["0-10", "10-20", "20-40", "40-80", "80-160", ">160"]
    prox["distance_bin"] = pd.cut(dist, bins=bins, labels=labels, right=False, include_lowest=True)
    prox.to_csv(OUT / "panel_K_observed_ec_distance_to_nearest_Met_hi_Fib.csv", index=False)
    frac = (
        prox.groupby("distance_bin", observed=True)["metEClabel"]
        .apply(lambda s: float((s.astype(str) == "Met_hi_EC").mean()))
        .rename("fraction_Met_hi_EC")
        .reset_index()
    )
    counts = prox.groupby("distance_bin", observed=True)["cell_id"].size().rename("n_ec").reset_index()
    frac = frac.merge(counts, on="distance_bin", how="left")
    frac.to_csv(OUT / "panel_K_fraction_Met_hi_EC_by_distance_bin.csv", index=False)

    set_theme("whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
    sns.boxplot(
        data=prox,
        x="metEClabel",
        y="distance_to_nearest_Met_hi_Fib",
        hue="metEClabel",
        palette=LABEL_COLORS,
        fliersize=0,
        width=0.55,
        ax=axes[0],
        legend=False,
    )
    sns.stripplot(
        data=prox,
        x="metEClabel",
        y="distance_to_nearest_Met_hi_Fib",
        color="black",
        alpha=0.35,
        size=2.8,
        jitter=0.18,
        ax=axes[0],
    )
    axes[0].set_xlabel("IMC EC label")
    axes[0].set_ylabel("Distance to nearest Met_hi_Fib")
    axes[0].set_title("Observed EC proximity to Met_hi_Fib")
    sns.barplot(data=frac, x="distance_bin", y="fraction_Met_hi_EC", color="#2b7bba", ax=axes[1])
    for i, row in frac.iterrows():
        axes[1].text(i, row["fraction_Met_hi_EC"] + 0.03, f"n={int(row['n_ec'])}", ha="center", fontsize=8)
    axes[1].set_ylim(0, 1.08)
    axes[1].set_xlabel("Distance bin to nearest Met_hi_Fib")
    axes[1].set_ylabel("Fraction Met_hi_EC")
    axes[1].tick_params(axis="x", rotation=30)
    axes[1].set_title("Met_hi_EC fraction by proximity")
    fig.tight_layout()
    return savefig(fig, "panel_K_observed_met_hi_ec_proximity_to_met_hi_fib")


def write_html_index(plot_paths):
    rows = []
    for p in sorted(OUT.glob("panel_*.png")):
        rows.append(f"<h3>{html.escape(p.stem)}</h3><img src='{html.escape(p.name)}' style='max-width:100%;'>")
    tables = "".join(
        f"<li>{html.escape(p.name)}</li>"
        for p in sorted(OUT.glob("panel_*.csv"))
    )
    content = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>Cross-modal validation of IMC-defined EC metabolic states</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 32px; color: #222; line-height: 1.45; }}
    img {{ border: 1px solid #ddd; margin-bottom: 24px; }}
    code {{ background: #f4f4f4; padding: 1px 4px; }}
  </style>
</head>
<body>
  <h1>Cross-modal validation of IMC-defined EC metabolic states</h1>
  <p>Input observed Xenium AnnData: <code>{html.escape(str(H5AD))}</code></p>
  <p>Output folder: <code>{html.escape(str(OUT))}</code></p>
  <h2>Tables</h2>
  <ul>{tables}</ul>
  <h2>Plots</h2>
  {''.join(rows)}
</body>
</html>
"""
    path = OUT / "cross_modal_ec_validation_report.html"
    path.write_text(content)
    return path


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    adata, ec, X_log, X_log_dense = load_data()
    paths = []
    paths.extend(panel_a_schematic())
    paths.extend(panel_b_embedding(ec, X_log_dense))
    meta, pb = pseudobulk_profiles(ec)
    meta.to_csv(OUT / "panel_C_pseudobulk_metadata.csv", index=False)
    pb.to_csv(OUT / "panel_C_pseudobulk_logcpm_matrix.csv", index=False)
    paths.extend(panel_c_pseudobulk_pca(meta, pb))
    de = differential_expression(ec, meta, pb)
    paths.extend(panel_d_volcano(de))
    paths.extend(panel_e_top_gene_heatmap_dotplot(ec, X_log_dense, de))
    score_df = compute_scores(ec, X_log_dense)
    paths.extend(panel_f_scores(score_df))
    paths.extend(panel_g_sample_paired_scores(score_df))
    paths.extend(panel_h_loso_classifier(ec, X_log))
    paths.extend(panel_i_observed_vs_generated_classifier(ec, X_log, adata))
    paths.extend(panel_j_spatial_scores(score_df, ec, X_log))
    paths.extend(panel_k_proximity(adata, ec))
    index = write_html_index(paths)
    for path in paths:
        print(path)
    print(index)


if __name__ == "__main__":
    main()
