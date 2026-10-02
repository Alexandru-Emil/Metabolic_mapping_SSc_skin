
import os
from pathlib import Path
def project_path(relative):
    return Path(os.environ.get("METABOLIC_INPUT_DIR", "inputs")) / relative

#!/usr/bin/env python

from pathlib import Path
import html
import json

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


BASE = Path(project_path('external/perturb_met_hi_fib_ec'))
OUT = Path(project_path('external/perturb_met_hi_fib_ec_imc_scores'))
CELL_RESULTS = BASE / "ec_cell_level_perturbation_results.csv"
IMC_SCORES = OUT / "imc_protein_scores_mapped_to_sfe_cells.csv"

SCENARIOS = ["baseline", "remove_Met_hi_Fib", "replace_Met_hi_Fib_with_Other_Fib"]
SCENARIO_LABELS = {
    "baseline": "Baseline",
    "remove_Met_hi_Fib": "Remove Met_hi_Fib",
    "replace_Met_hi_Fib_with_Other_Fib": "Replace Met_hi_Fib with Other_Fib",
}
STATE_ORDER = ["Met_hi_EC", "Met_int_EC"]
SCORE_COLS = [
    "imc_glycolysis_score",
    "imc_tca_score",
    "imc_hif1a_score",
    "imc_nox4_score",
]
SCORE_LABELS = {
    "imc_glycolysis_score": "IMC glycolysis",
    "imc_tca_score": "IMC TCA/OXPHOS",
    "imc_hif1a_score": "IMC HIF1a",
    "imc_nox4_score": "IMC NOX4",
}


def write_table(df, path):
    df.to_csv(path, index=False)


def summarize_frequencies(df):
    counts = (
        df.groupby(["scenario", "predicted_ec_state"], observed=True)["cell_id"]
        .size()
        .rename("n_ec_state")
        .reset_index()
    )
    totals = df.groupby("scenario", observed=True)["cell_id"].size().rename("n_ec_total").reset_index()
    full = pd.MultiIndex.from_product([SCENARIOS, STATE_ORDER], names=["scenario", "predicted_ec_state"])
    counts = counts.set_index(["scenario", "predicted_ec_state"]).reindex(full, fill_value=0).reset_index()
    out = counts.merge(totals, on="scenario", how="left")
    out["freq_ec_state"] = out["n_ec_state"] / out["n_ec_total"]
    out["scenario_label"] = out["scenario"].map(SCENARIO_LABELS)
    base = out[out["scenario"] == "baseline"][["predicted_ec_state", "freq_ec_state"]].rename(
        columns={"freq_ec_state": "baseline_freq_ec_state"}
    )
    out = out.merge(base, on="predicted_ec_state", how="left")
    out["delta_vs_baseline"] = out["freq_ec_state"] - out["baseline_freq_ec_state"]
    return out


def summarize_scores(df):
    long = df.melt(
        id_vars=["scenario", "sample_id", "cell_id", "predicted_ec_state", "metEClabel"],
        value_vars=SCORE_COLS,
        var_name="score",
        value_name="value",
    )
    summary = (
        long.groupby(["scenario", "predicted_ec_state", "score"], observed=True)
        .agg(n_cells=("cell_id", "size"), mean=("value", "mean"), median=("value", "median"))
        .reset_index()
    )
    sample = (
        long.groupby(["scenario", "sample_id", "predicted_ec_state", "score"], observed=True)
        .agg(n_cells=("cell_id", "size"), mean=("value", "mean"), median=("value", "median"))
        .reset_index()
    )
    return long, summary, sample


def plot_frequency(freq):
    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(10, 5.5))
    sns.barplot(data=freq, x="predicted_ec_state", y="freq_ec_state", hue="scenario_label", ax=ax)
    for container in ax.containers:
        ax.bar_label(container, fmt="%.2f", fontsize=9, padding=3)
    ax.set_ylim(0, 1.08)
    ax.set_xlabel("Predicted EC state")
    ax.set_ylabel("Frequency among all ECs")
    ax.set_title("MintFlow EC-state frequencies")
    ax.legend(title=None, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3)
    fig.tight_layout()
    fig.savefig(OUT / "imc_report_ec_state_frequencies.pdf")
    fig.savefig(OUT / "imc_report_ec_state_frequencies.png", dpi=220)
    plt.close(fig)


def plot_scores_by_predicted_state(long):
    sns.set_theme(style="whitegrid", context="talk")
    plot_df = long.copy()
    plot_df["scenario_label"] = plot_df["scenario"].map(SCENARIO_LABELS)
    plot_df["score_label"] = plot_df["score"].map(SCORE_LABELS)
    g = sns.catplot(
        data=plot_df,
        x="predicted_ec_state",
        y="value",
        hue="scenario_label",
        col="score_label",
        kind="bar",
        errorbar=("ci", 95),
        col_wrap=2,
        height=4.2,
        aspect=1.25,
        sharey=False,
    )
    g.set_axis_labels("Predicted EC state", "Mean IMC protein score")
    g.set_titles("{col_name}")
    sns.move_legend(g, "lower center", bbox_to_anchor=(0.5, -0.02), ncol=3, title=None)
    g.fig.subplots_adjust(bottom=0.16, top=0.9)
    g.fig.suptitle("IMC protein scores grouped by MintFlow-predicted EC state")
    g.fig.savefig(OUT / "imc_report_scores_by_predicted_ec_state.pdf", bbox_inches="tight")
    g.fig.savefig(OUT / "imc_report_scores_by_predicted_ec_state.png", dpi=220, bbox_inches="tight")
    plt.close(g.fig)


def plot_sample_scores(sample):
    sns.set_theme(style="whitegrid", context="talk")
    plot_df = sample.copy()
    plot_df["scenario_label"] = plot_df["scenario"].map(SCENARIO_LABELS)
    plot_df["score_label"] = plot_df["score"].map(SCORE_LABELS)
    g = sns.catplot(
        data=plot_df,
        x="scenario_label",
        y="mean",
        hue="predicted_ec_state",
        col="score_label",
        kind="strip",
        jitter=0.18,
        dodge=True,
        s=7,
        alpha=0.85,
        col_wrap=2,
        height=4.2,
        aspect=1.35,
        sharey=False,
    )
    for ax in g.axes.flat:
        ax.tick_params(axis="x", rotation=25)
        ax.set_xlabel("")
        ax.set_ylabel("Sample mean IMC protein score")
    g.set_titles("{col_name}")
    sns.move_legend(g, "lower center", bbox_to_anchor=(0.5, -0.02), ncol=2, title="Predicted EC state")
    g.fig.subplots_adjust(bottom=0.18, top=0.9)
    g.fig.suptitle("Sample-level IMC protein scores by predicted EC state")
    g.fig.savefig(OUT / "imc_report_sample_scores_by_predicted_ec_state.pdf", bbox_inches="tight")
    g.fig.savefig(OUT / "imc_report_sample_scores_by_predicted_ec_state.png", dpi=220, bbox_inches="tight")
    plt.close(g.fig)


def df_to_html(df, digits=4):
    return df.to_html(index=False, float_format=lambda x: f"{x:.{digits}g}")


def write_html(freq, score_summary, sample_summary):
    score_genes = {
        "glycolysis": ["GLUT1", "HK1", "PFKL1_PFKM", "PKM2", "LDHA", "GAPDH"],
        "tca_oxphos": ["CS", "OGDH", "ATP5A", "SDHA"],
        "hif1a": ["HIF1a"],
        "nox4": ["NOX4"],
        "assay": "IMC SPE assay zrescaled",
    }
    imgs = [
        "imc_report_ec_state_frequencies.png",
        "imc_report_scores_by_predicted_ec_state.png",
        "imc_report_sample_scores_by_predicted_ec_state.png",
    ]
    img_html = "".join(f"<h3>{html.escape(p)}</h3><img src='{html.escape(p)}'>" for p in imgs)
    page = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>MintFlow perturbation with IMC protein scores</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 32px; color: #222; line-height: 1.45; }}
    table {{ border-collapse: collapse; margin: 16px 0; font-size: 13px; }}
    th, td {{ border: 1px solid #ccc; padding: 5px 7px; text-align: right; }}
    th:first-child, td:first-child {{ text-align: left; }}
    img {{ max-width: 100%; height: auto; display: block; margin: 12px 0 28px 0; }}
    code {{ background: #f4f4f4; padding: 1px 4px; }}
    .note {{ background: #f7f7f7; border-left: 4px solid #777; padding: 12px 16px; }}
  </style>
</head>
<body>
  <h1>MintFlow Met_hi_Fib perturbation with IMC protein metabolic scores</h1>
  <div class="note">
    <p>The perturbation predictions are the existing Python MintFlow outputs from <code>{html.escape(str(BASE))}</code>.</p>
    <p>The metabolic scores in this report are <b>not</b> generated Xenium scores. They are observed IMC protein scores mapped to SFE/Xenium cells by <code>KEY</code>.</p>
    <p>Because IMC protein scores are observed baseline measurements, perturbation does not alter the protein values themselves. Perturbation changes how cells are grouped by MintFlow-predicted EC state.</p>
  </div>
  <h2>Protein score definitions</h2>
  <pre>{html.escape(json.dumps(score_genes, indent=2))}</pre>
  <h2>EC-state frequency summary</h2>
  {df_to_html(freq)}
  <h2>IMC score summary by predicted EC state</h2>
  {df_to_html(score_summary)}
  <h2>Sample-level IMC score summary</h2>
  {df_to_html(sample_summary)}
  <h2>Plots</h2>
  {img_html}
  <h2>Exported files</h2>
  <ul>
    <li>imc_score_cell_level_perturbation_results.csv</li>
    <li>imc_score_ec_state_frequency_summary.csv</li>
    <li>imc_score_summary_by_predicted_ec_state.csv</li>
    <li>imc_score_sample_summary_by_predicted_ec_state.csv</li>
  </ul>
</body>
</html>
"""
    (OUT / "mintflow_met_hi_fib_perturbation_ec_report_IMC_scores.html").write_text(page)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cell = pd.read_csv(CELL_RESULTS)
    scores = pd.read_csv(IMC_SCORES)
    df = cell.merge(scores, on="cell_id", how="left", suffixes=("", "_scoremap"))
    missing = df[SCORE_COLS].isna().any(axis=1).sum()
    if missing:
        raise RuntimeError(f"{missing} perturbation EC rows are missing IMC protein scores.")

    df["scenario"] = pd.Categorical(df["scenario"], SCENARIOS, ordered=True)
    df["predicted_ec_state"] = pd.Categorical(df["predicted_ec_state"], STATE_ORDER, ordered=True)

    freq = summarize_frequencies(df)
    long, score_summary, sample_summary = summarize_scores(df)

    write_table(df, OUT / "imc_score_cell_level_perturbation_results.csv")
    write_table(freq, OUT / "imc_score_ec_state_frequency_summary.csv")
    write_table(score_summary, OUT / "imc_score_summary_by_predicted_ec_state.csv")
    write_table(sample_summary, OUT / "imc_score_sample_summary_by_predicted_ec_state.csv")

    plot_frequency(freq)
    plot_scores_by_predicted_state(long)
    plot_sample_scores(sample_summary)
    write_html(freq, score_summary, sample_summary)
    print(OUT / "mintflow_met_hi_fib_perturbation_ec_report_IMC_scores.html")


if __name__ == "__main__":
    main()
