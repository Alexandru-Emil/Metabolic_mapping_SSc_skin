
import os
from pathlib import Path
def project_path(relative):
    return Path(os.environ.get("METABOLIC_INPUT_DIR", "inputs")) / relative

#!/usr/bin/env python3

from pathlib import Path
import html

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


OUT = Path(
    project_path('Revision/Nan annotation/Data/RMD files/MintFlow/07_main_mintflow_manuscript_figure')
)
SRC_CANDIDATES = [
    Path(project_path('external/perturb_met_hi_fib_ec_with_frequency_dotplots')),
    Path(project_path('Revision/Nan annotation/Data/RMD files/MintFlow/05_final_perturbation_report_and_paper_plots/perturb_met_hi_fib_ec_with_frequency_dotplots')),
]
H5AD_CANDIDATES = [
    Path(project_path('external/xenium_sfe_09032026_combined_fib_ec_metlabel.h5ad')),
    Path(project_path('Revision/Nan annotation/Data/RMD files/MintFlow/01_anndata_exports/xenium_sfe_09032026_combined_fib_ec_metlabel.h5ad')),
]

BASELINE = "baseline"
REPLACE = "replace_Met_hi_Fib_with_Other_Fib"
STATE_ORDER = ["Met_hi_EC-like", "Met_int_EC-like"]
LABEL_ORDER = ["Met_hi_EC", "Met_int_EC"]
STATE_PALETTE = {"Met_hi_EC-like": "#2b7bba", "Met_int_EC-like": "#55a868"}
LABEL_PALETTE = {"Met_hi_EC": "#2b7bba", "Met_int_EC": "#55a868"}

PAPER_RC = {
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 7,
    "legend.title_fontsize": 8,
}


def first_existing(paths, must_contain=None):
    for path in paths:
        if path.exists() and (must_contain is None or (path / must_contain).exists()):
            return path
    raise FileNotFoundError("None of these paths are available:\n" + "\n".join(map(str, paths)))


def set_theme(style="whitegrid"):
    sns.set_theme(style=style, context="paper", font_scale=1.15, rc=PAPER_RC)


def savefig(fig, stem, dpi=350):
    pdf = OUT / f"{stem}.pdf"
    png = OUT / f"{stem}.png"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return pdf, png


def load_data():
    src = first_existing(SRC_CANDIDATES, "ec_cell_level_perturbation_results.csv")
    h5ad = first_existing(H5AD_CANDIDATES)
    cell = pd.read_csv(src / "ec_cell_level_perturbation_results.csv")
    cell = cell[cell["scenario"].isin([BASELINE, REPLACE])].copy()
    cell["scenario_label"] = cell["scenario"].map(
        {BASELINE: "Baseline", REPLACE: "Replace Met_hi_Fib\nwith Other_Fib"}
    )
    cell["predicted_ec_state_like"] = cell["predicted_ec_state"].map(
        {"Met_hi_EC": "Met_hi_EC-like", "Met_int_EC": "Met_int_EC-like"}
    )

    paired_path = src / "paper_replacement_ec_paired_probability_delta.csv"
    if paired_path.exists():
        paired = pd.read_csv(paired_path)
    else:
        base = cell[cell["scenario"] == BASELINE][["cell_id", "metEClabel", "prob_Met_hi_EC"]].rename(
            columns={"prob_Met_hi_EC": "prob_baseline"}
        )
        repl = cell[cell["scenario"] == REPLACE][["cell_id", "prob_Met_hi_EC"]].rename(
            columns={"prob_Met_hi_EC": "prob_replace"}
        )
        paired = base.merge(repl, on="cell_id", how="inner")
        paired["delta_prob"] = paired["prob_replace"] - paired["prob_baseline"]

    adata = ad.read_h5ad(h5ad)
    adata.obs["cell_id"] = adata.obs["cell_id"].astype(str)
    adata.obs["combined_met_label"] = adata.obs["combined_met_label"].astype(str)
    coords = adata.obs[["cell_id", "x_centroid", "y_centroid", "combined_met_label"]].reset_index(drop=True)
    paired = paired.merge(coords, on="cell_id", how="left", suffixes=("", "_obs"))
    if "x_centroid_obs" in paired.columns:
        paired["x_centroid"] = paired["x_centroid"].fillna(paired["x_centroid_obs"])
        paired["y_centroid"] = paired["y_centroid"].fillna(paired["y_centroid_obs"])
    paired.to_csv(OUT / "main_figure_paired_ec_probability_delta.csv", index=False)
    return src, h5ad, cell, paired, adata


def draw_panel_a(ax):
    ax.axis("off")
    ax.set_title("A  Perturbation design", loc="left", fontweight="bold", pad=8)
    boxes = [
        (0.05, 0.75, 0.90, 0.15, "IMC-defined Met_hi_Fib\nprotein-derived label"),
        (0.05, 0.52, 0.90, 0.15, "Matched Xenium / MintFlow model\nRNA generative space"),
        (0.05, 0.29, 0.90, 0.15, "Replace Met_hi_Fib\nwith Other_Fib"),
        (0.05, 0.06, 0.90, 0.15, "Predict generated EC response\nMet_hi_EC-like probability"),
    ]
    for x, y, w, h, text in boxes:
        ax.add_patch(
            plt.Rectangle(
                (x, y),
                w,
                h,
                facecolor="#f7f7f7",
                edgecolor="#333333",
                linewidth=1.1,
                transform=ax.transAxes,
            )
        )
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=8.5, transform=ax.transAxes)
    for y0, y1 in [(0.75, 0.67), (0.52, 0.44), (0.29, 0.21)]:
        ax.annotate(
            "",
            xy=(0.5, y1),
            xytext=(0.5, y0),
            xycoords="axes fraction",
            arrowprops=dict(arrowstyle="->", linewidth=1.4, color="#333333"),
        )


def draw_panel_b(ax, cell):
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
        freq["scenario_label"], ["Baseline", "Replace Met_hi_Fib\nwith Other_Fib"], ordered=True
    )
    freq["predicted_ec_state_like"] = pd.Categorical(freq["predicted_ec_state_like"], STATE_ORDER, ordered=True)
    freq.to_csv(OUT / "main_figure_replacement_ec_state_frequency.csv", index=False)

    sns.barplot(
        data=freq,
        x="scenario_label",
        y="frequency",
        hue="predicted_ec_state_like",
        palette=STATE_PALETTE,
        ax=ax,
    )
    for container in ax.containers:
        ax.bar_label(container, labels=[f"{bar.get_height() * 100:.1f}%" for bar in container], fontsize=7, padding=2)
    ax.set_title("B  Replacement of Met_hi_Fib reduces Met_hi_EC-like generated ECs", loc="left", fontweight="bold", pad=8)
    ax.set_xlabel("")
    ax.set_ylabel("Predicted EC-state frequency")
    ax.set_ylim(0, 0.74)
    ax.legend(title="Generated EC identity", frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0))


def draw_panel_c(ax, paired):
    plot_df = paired.copy()
    plot_df["metEClabel"] = pd.Categorical(plot_df["metEClabel"], LABEL_ORDER, ordered=True)
    sns.boxplot(
        data=plot_df,
        x="metEClabel",
        y="delta_prob",
        hue="metEClabel",
        palette=LABEL_PALETTE,
        fliersize=0,
        width=0.55,
        ax=ax,
        legend=False,
    )
    sns.stripplot(
        data=plot_df,
        x="metEClabel",
        y="delta_prob",
        color="black",
        alpha=0.42,
        size=2.7,
        jitter=0.18,
        ax=ax,
    )
    ax.axhline(0, color="black", linestyle="--", linewidth=0.9)
    ax.set_title("C  Per-cell EC probability shift", loc="left", fontweight="bold", pad=8)
    ax.set_xlabel("Original IMC EC label")
    ax.set_ylabel("Delta probability\nreplacement - baseline")
    ax.set_ylim(-1.08, 1.08)


def draw_panel_d(fig, gs, paired, adata):
    axes = [fig.add_subplot(gs[0]), fig.add_subplot(gs[1])]
    fib = adata.obs[adata.obs["combined_met_label"] == "Met_hi_Fib"][["x_centroid", "y_centroid"]]
    vmax = float(np.nanquantile(np.abs(paired["delta_prob"]), 0.99))
    vmax = max(vmax, 0.1)
    sca = None
    for ax, label in zip(axes, LABEL_ORDER):
        cur = paired[paired["metEClabel"] == label].copy()
        ax.scatter(
            fib["x_centroid"],
            fib["y_centroid"],
            s=5,
            facecolors="none",
            edgecolors="#555555",
            linewidths=0.25,
            alpha=0.35,
            label="baseline Met_hi_Fib",
            rasterized=True,
        )
        sca = ax.scatter(
            cur["x_centroid"],
            cur["y_centroid"],
            c=cur["delta_prob"],
            cmap="coolwarm",
            vmin=-vmax,
            vmax=vmax,
            s=22,
            linewidths=0.15,
            edgecolors="black",
            alpha=0.96,
            rasterized=True,
        )
        ax.set_title(f"Original {label}", fontsize=9)
        ax.set_aspect("equal", adjustable="box")
        ax.invert_yaxis()
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    axes[0].text(
        0,
        1.08,
        "D  Spatial map of delta Met_hi_EC-like probability",
        transform=axes[0].transAxes,
        ha="left",
        va="bottom",
        fontsize=10,
        fontweight="bold",
    )
    cbar = fig.colorbar(sca, ax=axes, fraction=0.028, pad=0.02)
    cbar.set_label("Delta probability\nreplacement - baseline", fontsize=8)
    axes[1].legend(frameon=False, loc="lower right", fontsize=7)


def individual_panels(cell, paired, adata):
    set_theme("white")
    fig, ax = plt.subplots(figsize=(4.8, 4.0))
    draw_panel_a(ax)
    savefig(fig, "panel_A_perturbation_schematic")

    set_theme("whitegrid")
    fig, ax = plt.subplots(figsize=(6.7, 4.4))
    draw_panel_b(ax, cell)
    savefig(fig, "panel_B_replacement_ec_state_frequency")

    set_theme("whitegrid")
    fig, ax = plt.subplots(figsize=(4.6, 4.3))
    draw_panel_c(ax, paired)
    savefig(fig, "panel_C_per_cell_probability_shift_by_original_ec_state")

    set_theme("white")
    fig = plt.figure(figsize=(11.5, 5.2), constrained_layout=True)
    gs = fig.add_gridspec(1, 2)
    draw_panel_d(fig, [gs[0, 0], gs[0, 1]], paired, adata)
    savefig(fig, "panel_D_spatial_delta_map_faceted_by_original_ec_state")


def combined_figure(cell, paired, adata):
    set_theme("whitegrid")
    fig = plt.figure(figsize=(13.5, 10.5), constrained_layout=True)
    outer = fig.add_gridspec(2, 2, width_ratios=[0.9, 1.25], height_ratios=[0.9, 1.1])
    ax_a = fig.add_subplot(outer[0, 0])
    ax_b = fig.add_subplot(outer[0, 1])
    ax_c = fig.add_subplot(outer[1, 0])
    draw_panel_a(ax_a)
    draw_panel_b(ax_b, cell)
    draw_panel_c(ax_c, paired)
    d_grid = outer[1, 1].subgridspec(1, 2)
    draw_panel_d(fig, [d_grid[0, 0], d_grid[0, 1]], paired, adata)
    savefig(fig, "main_mintflow_replacement_4panel_figure", dpi=400)


def write_html(src, h5ad):
    plots = [
        "main_mintflow_replacement_4panel_figure.png",
        "panel_A_perturbation_schematic.png",
        "panel_B_replacement_ec_state_frequency.png",
        "panel_C_per_cell_probability_shift_by_original_ec_state.png",
        "panel_D_spatial_delta_map_faceted_by_original_ec_state.png",
    ]
    rows = "".join(
        f"<h3>{html.escape(p)}</h3><img src='{html.escape(p)}' style='max-width:100%;'>"
        for p in plots
    )
    content = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>Main MintFlow Met_hi_Fib Replacement Figure</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 32px; color: #222; line-height: 1.45; }}
    img {{ border: 1px solid #ddd; margin-bottom: 24px; }}
    code {{ background: #f4f4f4; padding: 1px 4px; }}
  </style>
</head>
<body>
  <h1>Main MintFlow Met_hi_Fib Replacement Figure</h1>
  <p>This report reorganizes the key MintFlow result from the frequency-dotplot report into a manuscript-focused 4-panel layout.</p>
  <p><b>Source report/data folder:</b> <code>{html.escape(str(src))}</code></p>
  <p><b>Observed Xenium H5AD:</b> <code>{html.escape(str(h5ad))}</code></p>
  <h2>Panels</h2>
  <ul>
    <li><b>A.</b> Perturbation schematic linking IMC protein-derived labels, matched Xenium, MintFlow replacement, and generated EC response.</li>
    <li><b>B.</b> Replacement EC-state frequency plot.</li>
    <li><b>C.</b> Per-cell delta probability, stratified by original IMC EC state.</li>
    <li><b>D.</b> Spatial delta map faceted by original IMC EC state, with baseline Met_hi_Fib positions overlaid.</li>
  </ul>
  {rows}
</body>
</html>
"""
    path = OUT / "main_mintflow_replacement_4panel_report.html"
    path.write_text(content)
    return path


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    src, h5ad, cell, paired, adata = load_data()
    individual_panels(cell, paired, adata)
    combined_figure(cell, paired, adata)
    report = write_html(src, h5ad)
    print(report)
    print(OUT / "main_mintflow_replacement_4panel_figure.pdf")
    print(OUT / "panel_D_spatial_delta_map_faceted_by_original_ec_state.pdf")


if __name__ == "__main__":
    main()
