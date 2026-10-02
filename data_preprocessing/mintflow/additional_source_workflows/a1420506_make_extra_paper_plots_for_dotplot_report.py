
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
import seaborn as sns
from scipy import sparse
from scipy.spatial import cKDTree
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


WORKDIR = Path(project_path('external/mintflow_xenium_training'))
BASE = WORKDIR / "perturb_met_hi_fib_ec"
OUT = WORKDIR / "perturb_met_hi_fib_ec_with_frequency_dotplots"
H5AD = WORKDIR / "xenium_sfe_09032026_combined_fib_ec_metlabel.h5ad"
HTML = OUT / "mintflow_met_hi_fib_perturbation_ec_report_with_frequency_dotplots.html"

BASELINE = "baseline"
REPLACE = "replace_Met_hi_Fib_with_Other_Fib"
SCENARIO_LABELS = {
    BASELINE: "Baseline",
    REPLACE: "Replace Met_hi_Fib\nwith Other_Fib",
}
EC_LABELS = ["Met_hi_EC", "Met_int_EC"]
EC_LIKE_LABELS = {"Met_hi_EC": "Met_hi_EC-like", "Met_int_EC": "Met_int_EC-like"}
PAPER_RC = {
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 9,
    "legend.title_fontsize": 10,
}


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


def train_ec_classifier_from_generated(adata, Xmic):
    ec_mask = adata.obs["metEClabel"].isin(EC_LABELS).to_numpy()
    X_train = as_log1p_nonnegative(Xmic[ec_mask, :])
    y_train = adata.obs.loc[ec_mask, "metEClabel"].astype(str).to_numpy()
    clf = make_pipeline(
        StandardScaler(with_mean=False),
        LogisticRegression(max_iter=5000, class_weight="balanced", solver="liblinear"),
    )
    clf.fit(X_train, y_train)
    return clf


def savefig(fig, stem):
    pdf = OUT / f"{stem}.pdf"
    png = OUT / f"{stem}.png"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return pdf, png


def load_inputs():
    cell = pd.read_csv(OUT / "ec_cell_level_perturbation_results.csv")
    cell = cell[cell["scenario"].isin([BASELINE, REPLACE])].copy()
    cell["scenario_label"] = cell["scenario"].map(SCENARIO_LABELS)
    cell["predicted_ec_state_like"] = cell["predicted_ec_state"].map(EC_LIKE_LABELS)

    adata = ad.read_h5ad(H5AD)
    adata.obs["cell_id"] = adata.obs["cell_id"].astype(str)
    adata.obs["combined_met_label"] = adata.obs["combined_met_label"].astype(str)
    adata.obs["metEClabel"] = adata.obs["metEClabel"].astype(str)
    return cell, adata


def make_replacement_frequency_plot(cell):
    freq = (
        cell.groupby(["scenario", "scenario_label", "predicted_ec_state_like"], observed=True)["cell_id"]
        .size()
        .rename("n")
        .reset_index()
    )
    totals = cell.groupby("scenario", observed=True)["cell_id"].size().rename("total").reset_index()
    freq = freq.merge(totals, on="scenario", how="left")
    freq["frequency"] = freq["n"] / freq["total"]
    freq["scenario_label"] = pd.Categorical(
        freq["scenario_label"],
        [SCENARIO_LABELS[BASELINE], SCENARIO_LABELS[REPLACE]],
        ordered=True,
    )

    sns.set_theme(style="whitegrid", context="paper", font_scale=1.2, rc=PAPER_RC)
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    sns.barplot(
        data=freq,
        x="scenario_label",
        y="frequency",
        hue="predicted_ec_state_like",
        palette=["#2b7bba", "#55a868"],
        ax=ax,
    )
    for container in ax.containers:
        ax.bar_label(container, labels=[f"{v.get_height() * 100:.1f}%" for v in container], fontsize=10, padding=3)
    ax.set_xlabel("")
    ax.set_ylabel("Predicted EC-state frequency")
    ax.set_ylim(0, 0.72)
    ax.legend(title="Generated EC identity", frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    ax.set_title("Replacement shifts generated EC identity")
    fig.tight_layout()
    return savefig(fig, "paper_replacement_only_ec_state_frequency")


def make_original_state_stratified_plot(cell):
    frac = (
        cell.groupby(["scenario", "scenario_label", "metEClabel"], observed=True)["predicted_ec_state"]
        .apply(lambda s: float((s == "Met_hi_EC").mean()))
        .rename("fraction_pred_met_hi")
        .reset_index()
    )
    frac["scenario_label"] = pd.Categorical(
        frac["scenario_label"],
        [SCENARIO_LABELS[BASELINE], SCENARIO_LABELS[REPLACE]],
        ordered=True,
    )
    frac["metEClabel"] = pd.Categorical(frac["metEClabel"], EC_LABELS, ordered=True)

    sns.set_theme(style="whitegrid", context="paper", font_scale=1.2, rc=PAPER_RC)
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    sns.pointplot(
        data=frac,
        x="scenario_label",
        y="fraction_pred_met_hi",
        hue="metEClabel",
        palette=["#2b7bba", "#55a868"],
        dodge=0.18,
        markers="o",
        linestyles="-",
        linewidth=2.0,
        markersize=8,
        ax=ax,
    )
    for _, row in frac.iterrows():
        x = 0 if row["scenario"] == BASELINE else 1
        offset = -0.07 if row["metEClabel"] == "Met_hi_EC" else 0.07
        ax.text(
            x + offset,
            row["fraction_pred_met_hi"] + 0.035,
            f"{100 * row['fraction_pred_met_hi']:.1f}%",
            ha="center",
            va="bottom",
            fontsize=9,
        )
    ax.set_xlabel("")
    ax.set_ylabel("Fraction predicted Met_hi_EC-like")
    ax.set_ylim(0, 1.08)
    ax.legend(title="Original IMC EC label", frameon=False, loc="center left", bbox_to_anchor=(1.01, 0.5))
    ax.set_title("Replacement response by original EC state")
    fig.tight_layout()
    return savefig(fig, "paper_original_ec_state_fraction_pred_met_hi")


def paired_probability_table(cell, adata):
    base = cell[cell["scenario"] == BASELINE][["cell_id", "metEClabel", "prob_Met_hi_EC"]].rename(
        columns={"prob_Met_hi_EC": "prob_baseline"}
    )
    repl = cell[cell["scenario"] == REPLACE][["cell_id", "prob_Met_hi_EC"]].rename(
        columns={"prob_Met_hi_EC": "prob_replace"}
    )
    coords = adata.obs[["cell_id", "x_centroid", "y_centroid"]].reset_index(drop=True)
    paired = base.merge(repl, on="cell_id", how="inner").merge(coords, on="cell_id", how="left")
    paired["delta_prob"] = paired["prob_replace"] - paired["prob_baseline"]

    fib = adata.obs[adata.obs["combined_met_label"] == "Met_hi_Fib"][["x_centroid", "y_centroid"]].to_numpy(float)
    ec_xy = paired[["x_centroid", "y_centroid"]].to_numpy(float)
    paired["distance_to_nearest_baseline_Met_hi_Fib"] = cKDTree(fib).query(ec_xy, k=1)[0]
    paired.to_csv(OUT / "paper_replacement_ec_paired_probability_delta.csv", index=False)
    return paired


def make_paired_probability_shift_plots(paired):
    sns.set_theme(style="whitegrid", context="paper", font_scale=1.2, rc=PAPER_RC)

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.3))
    sns.scatterplot(
        data=paired,
        x="prob_baseline",
        y="prob_replace",
        hue="metEClabel",
        palette=["#2b7bba", "#55a868"],
        s=35,
        alpha=0.78,
        linewidth=0,
        ax=axes[0],
    )
    axes[0].plot([0, 1], [0, 1], color="black", linewidth=1, linestyle="--")
    axes[0].set_xlabel("Baseline probability Met_hi_EC-like")
    axes[0].set_ylabel("Replacement probability Met_hi_EC-like")
    axes[0].set_xlim(-0.03, 1.03)
    axes[0].set_ylim(-0.03, 1.03)
    axes[0].legend(title="Original IMC EC label", frameon=False, loc="lower right")

    sns.boxplot(
        data=paired,
        x="metEClabel",
        y="delta_prob",
        hue="metEClabel",
        palette=["#2b7bba", "#55a868"],
        fliersize=0,
        width=0.55,
        ax=axes[1],
        legend=False,
    )
    sns.stripplot(
        data=paired,
        x="metEClabel",
        y="delta_prob",
        color="black",
        alpha=0.45,
        size=3,
        jitter=0.22,
        ax=axes[1],
    )
    axes[1].axhline(0, color="black", linewidth=1, linestyle="--")
    axes[1].set_xlabel("Original IMC EC label")
    axes[1].set_ylabel("Delta probability\nreplacement - baseline")
    axes[1].set_title("Per-cell replacement shift")
    fig.tight_layout()
    return savefig(fig, "paper_per_cell_paired_prob_met_hi_shift")


def make_spatial_delta_map(paired, adata):
    sns.set_theme(style="white", context="paper", font_scale=1.2, rc=PAPER_RC)
    vmax = float(np.nanquantile(np.abs(paired["delta_prob"]), 0.99))
    vmax = max(vmax, 0.05)

    fig, ax = plt.subplots(figsize=(11, 9))
    ax.scatter(
        adata.obs["x_centroid"],
        adata.obs["y_centroid"],
        s=0.8,
        c="#dddddd",
        linewidths=0,
        alpha=0.22,
        rasterized=True,
    )
    sca = ax.scatter(
        paired["x_centroid"],
        paired["y_centroid"],
        s=10,
        c=paired["delta_prob"],
        cmap="coolwarm",
        vmin=-vmax,
        vmax=vmax,
        linewidths=0,
        alpha=0.96,
        rasterized=True,
    )
    cbar = fig.colorbar(sca, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("Delta probability Met_hi_EC-like\nreplacement - baseline")
    ax.set_aspect("equal", adjustable="box")
    ax.invert_yaxis()
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Spatial map of EC response to Met_hi_Fib replacement")
    fig.tight_layout()
    return savefig(fig, "paper_spatial_delta_prob_met_hi_replacement_minus_baseline")


def make_distance_response_plot(paired):
    bins = [0, 10, 20, 40, 80, 160, np.inf]
    labels = ["0-10", "10-20", "20-40", "40-80", "80-160", ">160"]
    paired = paired.copy()
    paired["distance_bin_um"] = pd.cut(
        paired["distance_to_nearest_baseline_Met_hi_Fib"],
        bins=bins,
        labels=labels,
        include_lowest=True,
        right=False,
    )
    binned = (
        paired.groupby(["distance_bin_um", "metEClabel"], observed=True)["delta_prob"]
        .agg(median="median", q25=lambda s: s.quantile(0.25), q75=lambda s: s.quantile(0.75), n="size")
        .reset_index()
    )
    binned.to_csv(OUT / "paper_replacement_distance_to_met_hi_fib_binned_delta.csv", index=False)

    sns.set_theme(style="whitegrid", context="paper", font_scale=1.2, rc=PAPER_RC)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    sns.scatterplot(
        data=paired,
        x="distance_to_nearest_baseline_Met_hi_Fib",
        y="delta_prob",
        hue="metEClabel",
        palette=["#2b7bba", "#55a868"],
        s=30,
        alpha=0.68,
        linewidth=0,
        ax=axes[0],
    )
    axes[0].axhline(0, color="black", linewidth=1, linestyle="--")
    axes[0].set_xlabel("Distance to nearest baseline Met_hi_Fib")
    axes[0].set_ylabel("Delta probability\nreplacement - baseline")
    axes[0].legend(title="Original IMC EC label", frameon=False)

    sns.pointplot(
        data=paired,
        x="distance_bin_um",
        y="delta_prob",
        hue="metEClabel",
        palette=["#2b7bba", "#55a868"],
        errorbar=("pi", 50),
        dodge=0.25,
        markers="o",
        linestyles="-",
        ax=axes[1],
    )
    axes[1].axhline(0, color="black", linewidth=1, linestyle="--")
    axes[1].set_xlabel("Distance bin to nearest baseline Met_hi_Fib")
    axes[1].set_ylabel("Delta probability\nreplacement - baseline")
    axes[1].tick_params(axis="x", rotation=30)
    axes[1].legend(title="Original IMC EC label", frameon=False)
    fig.tight_layout()
    return savefig(fig, "paper_distance_to_met_hi_fib_vs_ec_delta_prob")


def make_realisation_uncertainty_plot(adata):
    with (BASE / "generation_baseline.pkl").open("rb") as handle:
        baseline_result = pickle.load(handle)
    with (BASE / "generation_replace_Met_hi_Fib_with_Other_Fib.pkl").open("rb") as handle:
        replace_result = pickle.load(handle)

    clf = train_ec_classifier_from_generated(adata, avg_generated_xmic(baseline_result))
    ec_mask = adata.obs["metEClabel"].isin(EC_LABELS).to_numpy()

    rows = []
    for scenario, result in [(BASELINE, baseline_result), (REPLACE, replace_result)]:
        for idx, realisation in enumerate(result["list_generated_realisations_ie_expressions"], start=1):
            X_ec = as_log1p_nonnegative(realisation["MintFLow_Generated_Xmic"][ec_mask, :])
            pred = clf.predict(X_ec)
            rows.append(
                {
                    "realisation": idx,
                    "scenario": scenario,
                    "scenario_label": SCENARIO_LABELS[scenario],
                    "fraction_pred_Met_hi_EC_like": float((pred == "Met_hi_EC").mean()),
                }
            )
    real = pd.DataFrame(rows)
    real["scenario_label"] = pd.Categorical(
        real["scenario_label"],
        [SCENARIO_LABELS[BASELINE], SCENARIO_LABELS[REPLACE]],
        ordered=True,
    )
    real.to_csv(OUT / "paper_replacement_realisation_level_frequency.csv", index=False)

    sns.set_theme(style="whitegrid", context="paper", font_scale=1.2, rc=PAPER_RC)
    fig, ax = plt.subplots(figsize=(6.5, 5.2))
    for rid, group in real.groupby("realisation", observed=True):
        group = group.sort_values("scenario_label")
        ax.plot(
            group["scenario_label"],
            group["fraction_pred_Met_hi_EC_like"],
            color="#999999",
            linewidth=1.1,
            alpha=0.75,
            zorder=1,
        )
    sns.stripplot(
        data=real,
        x="scenario_label",
        y="fraction_pred_Met_hi_EC_like",
        hue="scenario_label",
        palette=["#2b7bba", "#c44e52"],
        size=8,
        jitter=0.07,
        edgecolor="white",
        linewidth=0.6,
        legend=False,
        ax=ax,
        zorder=2,
    )
    ax.set_xlabel("")
    ax.set_ylabel("Fraction predicted Met_hi_EC-like")
    ax.set_ylim(0, 0.75)
    ax.set_title("Replacement effect across MintFlow realizations")
    fig.tight_layout()
    return savefig(fig, "paper_realisation_level_replacement_frequency_shift")


def update_html(plot_files):
    if not HTML.exists():
        return
    html_text = HTML.read_text()
    section_title = "<h2>Additional paper-focused replacement plots</h2>"
    items = []
    for png in plot_files:
        items.append(
            f"<h3>{html.escape(png.name)}</h3>"
            f"<img src='{html.escape(png.name)}' style='max-width: 100%;'>"
        )
    section = section_title + "".join(items)
    if section_title in html_text:
        before = html_text.split(section_title)[0]
        after = html_text.split("<h2>Run Parameters</h2>", 1)[-1]
        html_text = before + section + "<h2>Run Parameters</h2>" + after
    else:
        html_text = html_text.replace("<h2>Run Parameters</h2>", section + "<h2>Run Parameters</h2>", 1)
    HTML.write_text(html_text)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cell, adata = load_inputs()
    plot_paths = []

    plot_paths.extend(make_replacement_frequency_plot(cell))
    plot_paths.extend(make_original_state_stratified_plot(cell))
    paired = paired_probability_table(cell, adata)
    plot_paths.extend(make_paired_probability_shift_plots(paired))
    plot_paths.extend(make_spatial_delta_map(paired, adata))
    plot_paths.extend(make_distance_response_plot(paired))
    plot_paths.extend(make_realisation_uncertainty_plot(adata))

    pngs = [p for p in plot_paths if p.suffix == ".png"]
    update_html(pngs)
    for path in plot_paths:
        print(path)
    print(HTML)


if __name__ == "__main__":
    main()
