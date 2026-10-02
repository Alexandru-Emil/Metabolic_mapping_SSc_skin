
import os
from pathlib import Path
def project_path(relative):
    return Path(os.environ.get("METABOLIC_INPUT_DIR", "inputs")) / relative

#!/usr/bin/env python

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


SRC = Path(project_path('external/perturb_met_hi_fib_ec'))
OUT = Path(project_path('external/perturb_met_hi_fib_ec_with_frequency_dotplots'))
H5AD = Path(project_path('external/xenium_sfe_09032026_combined_fib_ec_metlabel.h5ad'))

SCENARIO_ORDER = ["baseline", "remove_Met_hi_Fib", "replace_Met_hi_Fib_with_Other_Fib"]
SCENARIO_TITLES = {
    "baseline": "Baseline",
    "remove_Met_hi_Fib": "Remove Met_hi_Fib",
    "replace_Met_hi_Fib_with_Other_Fib": "Replace Met_hi_Fib with Other_Fib",
}
PALETTE = {
    "Met_hi_Fib": "#c03a2b",
    "Other_Fib": "#7a7a7a",
    "Met_hi_EC": "#1f77b4",
    "Met_int_EC": "#2ca02c",
}


def load_tissues():
    adata = ad.read_h5ad(H5AD)
    adata.obs["cell_id"] = adata.obs["cell_id"].astype(str)
    adata.obs["combined_met_label"] = adata.obs["combined_met_label"].astype(str)
    tissues = {
        "baseline": adata.copy(),
        "remove_Met_hi_Fib": adata[adata.obs["combined_met_label"] != "Met_hi_Fib"].copy(),
        "replace_Met_hi_Fib_with_Other_Fib": adata.copy(),
    }
    tissues["replace_Met_hi_Fib_with_Other_Fib"].obs["combined_met_label"] = tissues[
        "replace_Met_hi_Fib_with_Other_Fib"
    ].obs["combined_met_label"].replace({"Met_hi_Fib": "Other_Fib"})
    return tissues


def finish_axes(ax):
    ax.set_aspect("equal", adjustable="box")
    ax.invert_yaxis()
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.tick_params(labelsize=9)


def save_highres(fig, stem):
    pdf = OUT / f"{stem}_highres.pdf"
    png = OUT / f"{stem}_highres.png"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, dpi=600, bbox_inches="tight")
    plt.close(fig)
    return pdf, png


def plot_celltype_maps(tissues):
    sns.set_theme(style="white", context="paper")
    fig, axes = plt.subplots(1, 3, figsize=(48, 16), sharex=True, sharey=True, constrained_layout=True)
    for ax, scenario in zip(axes, SCENARIO_ORDER):
        tissue = tissues[scenario]
        plot_df = tissue.obs[["combined_met_label", "x_centroid", "y_centroid"]].copy()
        plot_df["combined_met_label"] = plot_df["combined_met_label"].astype(str)
        for label, group in plot_df.groupby("combined_met_label", observed=True):
            ax.scatter(
                group["x_centroid"],
                group["y_centroid"],
                s=1.0,
                c=PALETTE.get(label, "#333333"),
                label=label,
                linewidths=0,
                alpha=0.78,
                rasterized=False,
            )
        ax.set_title(SCENARIO_TITLES[scenario], fontsize=18)
        finish_axes(ax)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="center right", frameon=False, fontsize=15, markerscale=6)
    return save_highres(fig, "tissue_maps_celltype_perturbations")


def plot_all_ec_probability(cell_df, tissues):
    sns.set_theme(style="white", context="paper")
    fig, axes = plt.subplots(1, 3, figsize=(48, 16), sharex=True, sharey=True, constrained_layout=True)
    sca = None
    for ax, scenario in zip(axes, SCENARIO_ORDER):
        tissue = tissues[scenario]
        ec = cell_df[cell_df["scenario"] == scenario].merge(
            tissue.obs[["cell_id", "x_centroid", "y_centroid"]].reset_index(drop=True),
            on="cell_id",
            how="left",
        )
        ax.scatter(
            tissue.obs["x_centroid"],
            tissue.obs["y_centroid"],
            s=0.45,
            c="#dddddd",
            linewidths=0,
            alpha=0.28,
            rasterized=False,
        )
        sca = ax.scatter(
            ec["x_centroid"],
            ec["y_centroid"],
            s=4.0,
            c=ec["prob_Met_hi_EC"],
            cmap="viridis",
            vmin=0,
            vmax=1,
            linewidths=0,
            alpha=0.95,
            rasterized=False,
        )
        ax.set_title(SCENARIO_TITLES[scenario], fontsize=18)
        finish_axes(ax)
    cbar = fig.colorbar(sca, ax=axes, fraction=0.025, pad=0.015)
    cbar.set_label("Generated EC probability: Met_hi_EC", fontsize=15)
    cbar.ax.tick_params(labelsize=12)
    return save_highres(fig, "tissue_maps_ec_met_hi_probability")


def plot_original_ec_probability(cell_df, tissues, label):
    sns.set_theme(style="white", context="paper")
    fig, axes = plt.subplots(1, 3, figsize=(48, 16), sharex=True, sharey=True, constrained_layout=True)
    sca = None
    for ax, scenario in zip(axes, SCENARIO_ORDER):
        tissue = tissues[scenario]
        ec = cell_df[(cell_df["scenario"] == scenario) & (cell_df["metEClabel"] == label)].merge(
            tissue.obs[["cell_id", "x_centroid", "y_centroid"]].reset_index(drop=True),
            on="cell_id",
            how="left",
        )
        ax.scatter(
            tissue.obs["x_centroid"],
            tissue.obs["y_centroid"],
            s=0.45,
            c="#e5e5e5",
            linewidths=0,
            alpha=0.22,
            rasterized=False,
        )
        sca = ax.scatter(
            ec["x_centroid"],
            ec["y_centroid"],
            s=5.0,
            c=ec["prob_Met_hi_EC"],
            cmap="viridis",
            vmin=0,
            vmax=1,
            linewidths=0,
            alpha=0.98,
            rasterized=False,
        )
        ax.set_title(f"{SCENARIO_TITLES[scenario]}: original {label}", fontsize=18)
        finish_axes(ax)
    cbar = fig.colorbar(sca, ax=axes, fraction=0.025, pad=0.015)
    cbar.set_label("Generated EC probability: Met_hi_EC", fontsize=15)
    cbar.ax.tick_params(labelsize=12)
    return save_highres(fig, f"tissue_maps_{label.lower()}_met_hi_probability")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cell_df = pd.read_csv(OUT / "ec_cell_level_perturbation_results.csv")
    tissues = load_tissues()
    outputs = []
    outputs.extend(plot_celltype_maps(tissues))
    outputs.extend(plot_all_ec_probability(cell_df, tissues))
    outputs.extend(plot_original_ec_probability(cell_df, tissues, "Met_hi_EC"))
    outputs.extend(plot_original_ec_probability(cell_df, tissues, "Met_int_EC"))
    for path in outputs:
        print(path)


if __name__ == "__main__":
    main()
