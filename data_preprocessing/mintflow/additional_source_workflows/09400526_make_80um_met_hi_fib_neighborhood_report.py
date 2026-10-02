
import os
from pathlib import Path
def project_path(relative):
    return Path(os.environ.get("METABOLIC_INPUT_DIR", "inputs")) / relative

#!/usr/bin/env python

import argparse
from pathlib import Path
import re

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu


BASE_DIR = Path(project_path('external/perturb_met_hi_fib_ec'))
H5AD = Path(project_path('external/xenium_sfe_09032026_combined_fib_ec_metlabel.h5ad'))
CELL_RESULTS = BASE_DIR / "ec_cell_level_perturbation_results.csv"

SCENARIO_ORDER = ["baseline", "remove_Met_hi_Fib", "replace_Met_hi_Fib_with_Other_Fib"]
SCENARIO_LABELS = {
    "baseline": "Baseline",
    "remove_Met_hi_Fib": "Remove Met_hi_Fib",
    "replace_Met_hi_Fib_with_Other_Fib": "Replace Met_hi_Fib with Other_Fib",
}
STATE_ORDER = ["Met_hi_EC", "Met_int_EC"]


def ensure_dir(path: Path):
    path.mkdir(parents=True, exist_ok=True)


def radius_tag(radius_um: float) -> str:
    if float(radius_um).is_integer():
        return f"{int(radius_um)}um"
    return re.sub(r"[^0-9A-Za-z]+", "p", f"{radius_um:g}") + "um"


def compute_neighborhood(adata: ad.AnnData, radius_um: float) -> pd.DataFrame:
    obs = adata.obs.copy()
    obs["cell_id"] = obs["cell_id"].astype(str)
    ec = obs[obs["metEClabel"].isin(STATE_ORDER)].copy()
    fib = obs[obs["combined_met_label"] == "Met_hi_Fib"].copy()

    if fib.empty:
        raise RuntimeError("No baseline Met_hi_Fib cells found in the AnnData object.")
    if ec.empty:
        raise RuntimeError("No baseline EC cells found in the AnnData object.")

    fib_xy = fib[["x_centroid", "y_centroid"]].to_numpy(dtype=float)
    ec_xy = ec[["x_centroid", "y_centroid"]].to_numpy(dtype=float)

    tree = cKDTree(fib_xy)
    dist, idx = tree.query(ec_xy, k=1)
    ec["nearest_met_hi_fib_distance_um"] = dist
    ec["nearest_met_hi_fib_cell_id"] = fib.iloc[idx]["cell_id"].to_numpy()
    ec[f"within_{radius_tag(radius_um)}_met_hi_fib"] = ec["nearest_met_hi_fib_distance_um"] <= radius_um
    return ec[
        [
            "cell_id",
            "sample_id",
            "metEClabel",
            "x_centroid",
            "y_centroid",
            "nearest_met_hi_fib_distance_um",
            "nearest_met_hi_fib_cell_id",
            f"within_{radius_tag(radius_um)}_met_hi_fib",
        ]
    ].reset_index(drop=True)


def summarize_frequencies(cell_results: pd.DataFrame, neighborhood: pd.DataFrame, radius_um: float, denominator_mode: str):
    flag_col = f"within_{radius_tag(radius_um)}_met_hi_fib"
    subset_ids = set(neighborhood.loc[neighborhood[flag_col], "cell_id"])
    df = cell_results[cell_results["cell_id"].isin(subset_ids)].copy()
    df["scenario"] = pd.Categorical(df["scenario"], SCENARIO_ORDER, ordered=True)
    df["predicted_ec_state"] = pd.Categorical(df["predicted_ec_state"], STATE_ORDER, ordered=True)
    all_ec_results = cell_results.copy()
    all_ec_results["scenario"] = pd.Categorical(all_ec_results["scenario"], SCENARIO_ORDER, ordered=True)

    if denominator_mode == "all_ecs":
        totals = all_ec_results.groupby("scenario", observed=True)["cell_id"].size().rename("n_ec_total").reset_index()
    else:
        totals = df.groupby("scenario", observed=True)["cell_id"].size().rename("n_ec_total").reset_index()
    counts = (
        df.groupby(["scenario", "predicted_ec_state"], observed=True)["cell_id"]
        .size()
        .rename("n_ec_state")
        .reset_index()
    )
    full_index = pd.MultiIndex.from_product(
        [pd.Categorical(SCENARIO_ORDER, SCENARIO_ORDER, ordered=True), STATE_ORDER],
        names=["scenario", "predicted_ec_state"],
    )
    counts = counts.set_index(["scenario", "predicted_ec_state"]).reindex(full_index, fill_value=0).reset_index()
    summary = counts.merge(totals, on="scenario", how="left")
    summary["freq_ec_state"] = summary["n_ec_state"] / summary["n_ec_total"]
    summary["scenario_label"] = summary["scenario"].astype(str).map(SCENARIO_LABELS)

    if denominator_mode == "all_ecs":
        sample_totals = (
            all_ec_results.groupby(["scenario", "sample_id"], observed=True)["cell_id"]
            .size()
            .rename("n_ec_total")
            .reset_index()
        )
    else:
        sample_totals = df.groupby(["scenario", "sample_id"], observed=True)["cell_id"].size().rename("n_ec_total").reset_index()
    sample_counts = (
        df.groupby(["scenario", "sample_id", "predicted_ec_state"], observed=True)["cell_id"]
        .size()
        .rename("n_ec_state")
        .reset_index()
    )
    samples = sorted(all_ec_results["sample_id"].unique() if denominator_mode == "all_ecs" else df["sample_id"].unique())
    sample_index = pd.MultiIndex.from_product(
        [pd.Categorical(SCENARIO_ORDER, SCENARIO_ORDER, ordered=True), samples, STATE_ORDER],
        names=["scenario", "sample_id", "predicted_ec_state"],
    )
    sample_counts = (
        sample_counts.set_index(["scenario", "sample_id", "predicted_ec_state"])
        .reindex(sample_index, fill_value=0)
        .reset_index()
    )
    sample_summary = sample_counts.merge(sample_totals, on=["scenario", "sample_id"], how="left")
    sample_summary["freq_ec_state"] = sample_summary["n_ec_state"] / sample_summary["n_ec_total"]
    sample_summary["scenario_label"] = sample_summary["scenario"].astype(str).map(SCENARIO_LABELS)

    baseline = summary[summary["scenario"].astype(str) == "baseline"][
        ["predicted_ec_state", "freq_ec_state"]
    ].rename(columns={"freq_ec_state": "baseline_freq_ec_state"})
    summary = summary.merge(baseline, on="predicted_ec_state", how="left")
    summary["delta_vs_baseline"] = summary["freq_ec_state"] - summary["baseline_freq_ec_state"]
    return df, summary, sample_summary


def denominator_label(denominator_mode: str):
    return "all ECs" if denominator_mode == "all_ecs" else "local ECs"


def plot_overall(summary: pd.DataFrame, out_dir: Path, radius_um: float, denominator_mode: str):
    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(10, 5.5))
    sns.barplot(
        data=summary,
        x="predicted_ec_state",
        y="freq_ec_state",
        hue="scenario_label",
        ax=ax,
    )
    for container in ax.containers:
        ax.bar_label(container, fmt="%.2f", fontsize=9, padding=3)
    ax.set_ylim(0, 1.08)
    ax.set_xlabel("Predicted EC state")
    ax.set_ylabel(f"Local EC-state count / {denominator_label(denominator_mode)}")
    ax.set_title("EC-state frequencies near baseline Met_hi_Fib")
    ax.legend(title=None, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3)
    fig.tight_layout()
    tag = radius_tag(radius_um)
    suffix = "all_ec_denominator" if denominator_mode == "all_ecs" else "local_denominator"
    fig.savefig(out_dir / f"ec_state_frequency_within_{tag}_met_hi_fib_{suffix}_overall.pdf")
    fig.savefig(out_dir / f"ec_state_frequency_within_{tag}_met_hi_fib_{suffix}_overall.png", dpi=220)
    plt.close(fig)


def add_median_iqr(ax, data: pd.DataFrame, x_order):
    stats = data.groupby("scenario_label", observed=True)["freq_ec_state"].agg(
        median="median",
        q25=lambda s: s.quantile(0.25),
        q75=lambda s: s.quantile(0.75),
    )
    for i, label in enumerate(x_order):
        if label not in stats.index:
            continue
        row = stats.loc[label]
        ax.vlines(i, row["q25"], row["q75"], color="black", linewidth=2.2, zorder=4)
        ax.hlines(row["median"], i - 0.28, i + 0.28, color="black", linewidth=3, zorder=5)


def format_p(p):
    if pd.isna(p):
        return "MWU p = NA"
    if p < 0.001:
        return "MWU p < 0.001"
    return f"MWU p = {p:.3f}"


def plot_sample_pair(sample_summary: pd.DataFrame, out_dir: Path, radius_um: float, target_scenario: str, denominator_mode: str):
    sns.set_theme(style="whitegrid", context="talk")
    pair_order = ["baseline", target_scenario]
    x_order = [SCENARIO_LABELS[s] for s in pair_order]
    plot_df = sample_summary[sample_summary["scenario"].astype(str).isin(pair_order)].copy()
    plot_df["scenario_label"] = pd.Categorical(plot_df["scenario_label"], x_order, ordered=True)
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 6.2), sharey=True)
    for ax, state in zip(axes, STATE_ORDER):
        state_df = plot_df[plot_df["predicted_ec_state"].astype(str) == state]
        sns.stripplot(
            data=state_df,
            x="scenario_label",
            y="freq_ec_state",
            order=x_order,
            color="#334e68",
            size=11,
            jitter=0.24,
            alpha=0.82,
            edgecolor="white",
            linewidth=0.6,
            ax=ax,
        )
        add_median_iqr(ax, state_df, x_order)
        vals0 = state_df[state_df["scenario"].astype(str) == "baseline"]["freq_ec_state"]
        vals1 = state_df[state_df["scenario"].astype(str) == target_scenario]["freq_ec_state"]
        if len(vals0) and len(vals1):
            p = mannwhitneyu(vals0, vals1, alternative="two-sided").pvalue
        else:
            p = np.nan
        ax.text(
            0.5,
            1.12,
            format_p(p),
            ha="center",
            va="bottom",
            fontsize=11,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="#bbbbbb", alpha=0.95),
            clip_on=False,
        )
        ax.set_ylim(-0.03, 1.18)
        ax.tick_params(axis="x", rotation=25)
        ax.set_xlabel("")
        ax.set_ylabel(f"Local EC-state count / {denominator_label(denominator_mode)}\nwithin {radius_um:g} um" if ax is axes[0] else "")
        ax.set_title(state)
    fig.suptitle(f"Per-sample local EC-state frequencies: Baseline vs {SCENARIO_LABELS[target_scenario]}", y=1.04)
    fig.subplots_adjust(top=0.78, bottom=0.22, wspace=0.16)
    tag = radius_tag(radius_um)
    suffix = "all_ec_denominator" if denominator_mode == "all_ecs" else "local_denominator"
    short = "remove" if target_scenario == "remove_Met_hi_Fib" else "replace"
    fig.savefig(out_dir / f"ec_state_frequency_within_{tag}_met_hi_fib_{suffix}_per_sample_baseline_vs_{short}.pdf", bbox_inches="tight")
    fig.savefig(out_dir / f"ec_state_frequency_within_{tag}_met_hi_fib_{suffix}_per_sample_baseline_vs_{short}.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


def plot_sample(sample_summary: pd.DataFrame, out_dir: Path, radius_um: float, denominator_mode: str):
    plot_sample_pair(sample_summary, out_dir, radius_um, "remove_Met_hi_Fib", denominator_mode)
    plot_sample_pair(sample_summary, out_dir, radius_um, "replace_Met_hi_Fib_with_Other_Fib", denominator_mode)


def plot_distance_distribution(neighborhood: pd.DataFrame, out_dir: Path, radius_um: float):
    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.histplot(
        data=neighborhood,
        x="nearest_met_hi_fib_distance_um",
        hue="metEClabel",
        bins=30,
        element="step",
        stat="count",
        ax=ax,
    )
    ax.axvline(radius_um, color="black", linestyle="--", linewidth=1.5)
    ax.set_xlabel("Distance to nearest baseline Met_hi_Fib (um)")
    ax.set_ylabel("Number of ECs")
    ax.set_title("EC distance to nearest baseline Met_hi_Fib")
    fig.tight_layout()
    tag = radius_tag(radius_um)
    fig.savefig(out_dir / f"ec_distance_to_nearest_met_hi_fib_within_{tag}_cutoff.pdf")
    fig.savefig(out_dir / f"ec_distance_to_nearest_met_hi_fib_within_{tag}_cutoff.png", dpi=220)
    plt.close(fig)


def write_html(md_path: Path, html_path: Path):
    try:
        import markdown
        body = markdown.markdown(md_path.read_text(), extensions=["tables"])
    except Exception:
        body = "<pre>" + md_path.read_text() + "</pre>"
    html = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>{md_path.stem}</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 32px; color: #222; line-height: 1.45; }}
    table {{ border-collapse: collapse; margin: 16px 0; font-size: 13px; }}
    th, td {{ border: 1px solid #ccc; padding: 5px 7px; text-align: right; }}
    th:first-child, td:first-child {{ text-align: left; }}
    img {{ max-width: 100%; height: auto; display: block; margin: 12px 0 28px 0; }}
    code {{ background: #f4f4f4; padding: 1px 4px; }}
  </style>
</head>
<body>
{body}
</body>
</html>
"""
    html_path.write_text(html)


def write_markdown(neighborhood, summary, sample_summary, out_dir: Path, radius_um: float, denominator_mode: str):
    n_ec = len(neighborhood)
    tag = radius_tag(radius_um)
    suffix = "all_ec_denominator" if denominator_mode == "all_ecs" else "local_denominator"
    flag_col = f"within_{tag}_met_hi_fib"
    n_local = int(neighborhood[flag_col].sum())
    state_counts = (
        neighborhood[neighborhood[flag_col]]
        .groupby("metEClabel")
        .size()
        .reindex(STATE_ORDER, fill_value=0)
    )

    summary_md = summary.copy()
    summary_md["freq_ec_state"] = summary_md["freq_ec_state"].map(lambda x: f"{x:.4f}")
    summary_md["delta_vs_baseline"] = summary_md["delta_vs_baseline"].map(lambda x: f"{x:.4f}")

    sample_md = sample_summary.copy()
    sample_md["freq_ec_state"] = sample_md["freq_ec_state"].map(lambda x: f"{x:.4f}")

    def markdown_table(df):
        df = df.copy()
        df.columns = [str(c) for c in df.columns]
        lines = []
        lines.append("| " + " | ".join(df.columns) + " |")
        lines.append("| " + " | ".join(["---"] * len(df.columns)) + " |")
        for _, row in df.iterrows():
            lines.append("| " + " | ".join(str(row[c]) for c in df.columns) + " |")
        return "\n".join(lines)

    remove_plot = f"ec_state_frequency_within_{tag}_met_hi_fib_{suffix}_per_sample_baseline_vs_remove.png"
    replace_plot = f"ec_state_frequency_within_{tag}_met_hi_fib_{suffix}_per_sample_baseline_vs_replace.png"
    overall_plot = f"ec_state_frequency_within_{tag}_met_hi_fib_{suffix}_overall.png"
    distance_plot = f"ec_distance_to_nearest_met_hi_fib_within_{tag}_cutoff.png"

    md = f"""# Python MintFlow Perturbation: EC Frequencies Within {radius_um:g} um of Baseline Met_hi_Fib

This is a separate Python-based analysis using the existing MintFlow perturbation outputs. It does not use the R report pipeline.

## Definition

Local ECs were defined from the baseline tissue geometry as ECs whose centroid is within `{radius_um:g} um` of the nearest baseline `Met_hi_Fib` centroid.

- Total ECs in the trained AnnData: `{n_ec}`
- ECs within `{radius_um:g} um` of baseline `Met_hi_Fib`: `{n_local}`
- Baseline local EC annotations: `Met_hi_EC = {int(state_counts['Met_hi_EC'])}`, `Met_int_EC = {int(state_counts['Met_int_EC'])}`
- Frequency denominator: `{denominator_label(denominator_mode)}`

The same local EC cell IDs were then followed across baseline, `remove_Met_hi_Fib`, and `replace_Met_hi_Fib_with_Other_Fib` generated-expression predictions.

## Overall Frequencies

{markdown_table(summary_md[['scenario_label', 'predicted_ec_state', 'n_ec_state', 'n_ec_total', 'freq_ec_state', 'delta_vs_baseline']])}

![Overall local EC frequencies]({overall_plot})

## Per-Sample Frequencies

Each point is one sample. The numerator is the number of ECs within `{radius_um:g} um` of baseline `Met_hi_Fib` predicted as the given EC state; the denominator is `{denominator_label(denominator_mode)}` in that sample. Dots are deliberately not colored by sample. Black vertical lines show the IQR, black horizontal bars show the median, and each panel reports a two-sided Mann-Whitney U test.

### Baseline vs Remove Met_hi_Fib

![Per-sample local EC frequencies baseline vs remove]({remove_plot})

### Baseline vs Replace Met_hi_Fib with Other_Fib

![Per-sample local EC frequencies baseline vs replace]({replace_plot})

## Distance Check

![EC distance to nearest Met_hi_Fib]({distance_plot})

## Output Tables

- `ec_within_{tag}_met_hi_fib_cell_ids.csv`: EC IDs, sample, original EC label, and distance to nearest baseline `Met_hi_Fib`.
- `ec_within_{tag}_met_hi_fib_cell_level_predictions.csv`: scenario-level predictions for the local EC subset.
- `ec_within_{tag}_met_hi_fib_{suffix}_frequency_summary.csv`: overall frequency table.
- `ec_within_{tag}_met_hi_fib_{suffix}_sample_frequency_summary.csv`: per-sample frequency table.

## Plot Exports

- `ec_state_frequency_within_{tag}_met_hi_fib_{suffix}_overall.pdf`
- `ec_state_frequency_within_{tag}_met_hi_fib_{suffix}_per_sample_baseline_vs_remove.pdf`
- `ec_state_frequency_within_{tag}_met_hi_fib_{suffix}_per_sample_baseline_vs_replace.pdf`
- `ec_distance_to_nearest_met_hi_fib_within_{tag}_cutoff.pdf`
"""
    md_path = out_dir / f"ec_within_{tag}_met_hi_fib_{suffix}_report.md"
    html_path = out_dir / f"ec_within_{tag}_met_hi_fib_{suffix}_report.html"
    md_path.write_text(md)
    write_html(md_path, html_path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--radius-um", type=float, default=80.0)
    parser.add_argument("--denominator-mode", choices=["local_ecs", "all_ecs"], default="local_ecs")
    parser.add_argument("--outdir", type=Path, default=None)
    args = parser.parse_args()

    tag = radius_tag(args.radius_um)
    suffix = "all_ec_denominator" if args.denominator_mode == "all_ecs" else None
    out_dir = args.outdir or (BASE_DIR / (f"within_{tag}_met_hi_fib_{suffix}" if suffix else f"within_{tag}_met_hi_fib"))
    ensure_dir(out_dir)
    adata = ad.read_h5ad(H5AD)
    neighborhood = compute_neighborhood(adata, args.radius_um)
    cell_results = pd.read_csv(CELL_RESULTS)
    local_predictions, summary, sample_summary = summarize_frequencies(
        cell_results,
        neighborhood,
        args.radius_um,
        args.denominator_mode,
    )

    neighborhood.to_csv(out_dir / f"ec_within_{tag}_met_hi_fib_cell_ids.csv", index=False)
    local_predictions.to_csv(out_dir / f"ec_within_{tag}_met_hi_fib_cell_level_predictions.csv", index=False)
    file_suffix = "all_ec_denominator" if args.denominator_mode == "all_ecs" else "local_denominator"
    summary.to_csv(out_dir / f"ec_within_{tag}_met_hi_fib_{file_suffix}_frequency_summary.csv", index=False)
    sample_summary.to_csv(out_dir / f"ec_within_{tag}_met_hi_fib_{file_suffix}_sample_frequency_summary.csv", index=False)

    plot_overall(summary, out_dir, args.radius_um, args.denominator_mode)
    plot_sample(sample_summary, out_dir, args.radius_um, args.denominator_mode)
    plot_distance_distribution(neighborhood, out_dir, args.radius_um)
    write_markdown(neighborhood, summary, sample_summary, out_dir, args.radius_um, args.denominator_mode)
    print(f"Wrote neighborhood report to {out_dir}")
    print(summary)


if __name__ == "__main__":
    main()
