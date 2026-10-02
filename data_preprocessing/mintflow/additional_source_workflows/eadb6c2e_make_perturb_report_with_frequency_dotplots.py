
import os
from pathlib import Path
def project_path(relative):
    return Path(os.environ.get("METABOLIC_INPUT_DIR", "inputs")) / relative

#!/usr/bin/env python

from pathlib import Path
import re
import shutil

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


SRC = Path(project_path('external/perturb_met_hi_fib_ec'))
OUT = Path(project_path('external/perturb_met_hi_fib_ec_with_frequency_dotplots'))

SCENARIO_ORDER = ["baseline", "remove_Met_hi_Fib", "replace_Met_hi_Fib_with_Other_Fib"]
SCENARIO_LABELS = {
    "baseline": "Baseline",
    "remove_Met_hi_Fib": "Remove Met_hi_Fib",
    "replace_Met_hi_Fib_with_Other_Fib": "Replace Met_hi_Fib with Other_Fib",
}
STATE_ORDER = ["Met_hi_EC", "Met_int_EC"]


def copy_assets():
    OUT.mkdir(parents=True, exist_ok=True)
    patterns = ["*.csv", "*.json", "*.png", "*.pdf"]
    for pattern in patterns:
        for path in SRC.glob(pattern):
            if path.name.startswith("generation_"):
                continue
            shutil.copy2(path, OUT / path.name)


def make_pairwise_dotplot():
    freq = pd.read_csv(SRC / "ec_state_frequency_among_all_ecs.csv")
    freq = freq[freq["analysis_set"] == "all_ecs"].copy()
    freq["scenario_label"] = freq["scenario"].map(SCENARIO_LABELS)
    freq["scenario_label"] = pd.Categorical(
        freq["scenario_label"],
        [SCENARIO_LABELS[s] for s in SCENARIO_ORDER],
        ordered=True,
    )
    freq["predicted_ec_state"] = pd.Categorical(freq["predicted_ec_state"], STATE_ORDER, ordered=True)

    baseline = freq[freq["scenario"] == "baseline"].copy()
    remove = pd.concat([baseline, freq[freq["scenario"] == "remove_Met_hi_Fib"]], ignore_index=True)
    remove["comparison"] = "Baseline vs remove Met_hi_Fib"
    remove["comparison_scenario"] = remove["scenario"].map(
        {
            "baseline": "Baseline",
            "remove_Met_hi_Fib": "Remove",
        }
    )
    replace = pd.concat([baseline, freq[freq["scenario"] == "replace_Met_hi_Fib_with_Other_Fib"]], ignore_index=True)
    replace["comparison"] = "Baseline vs replace Met_hi_Fib"
    replace["comparison_scenario"] = replace["scenario"].map(
        {
            "baseline": "Baseline",
            "replace_Met_hi_Fib_with_Other_Fib": "Replace",
        }
    )
    plot_df = pd.concat([remove, replace], ignore_index=True)

    sns.set_theme(style="white", context="talk")
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2), sharey=True, constrained_layout=True)
    comparisons = ["Baseline vs remove Met_hi_Fib", "Baseline vs replace Met_hi_Fib"]
    cmap = sns.color_palette("viridis", as_cmap=True)
    norm = plt.Normalize(
        max(0, plot_df["freq_ec_state"].min() - 0.03),
        min(1, plot_df["freq_ec_state"].max() + 0.03),
    )
    for ax, comparison in zip(axes, comparisons):
        cur = plot_df[plot_df["comparison"] == comparison].copy()
        cur["x"] = cur["comparison_scenario"].map({"Baseline": 0, "Remove": 1, "Replace": 1})
        cur["y"] = cur["predicted_ec_state"].map({state: i for i, state in enumerate(STATE_ORDER)})
        ax.scatter(
            cur["x"],
            cur["y"],
            c=cur["freq_ec_state"],
            s=950,
            cmap=cmap,
            norm=norm,
            edgecolors="black",
            linewidths=0.8,
        )
        for _, row in cur.iterrows():
            ax.text(
                row["x"],
                row["y"],
                f"{100 * row['freq_ec_state']:.1f}%",
                ha="center",
                va="center",
                fontsize=11,
                color="white",
                fontweight="bold",
            )
        ax.set_title(comparison)
        ax.set_xlabel("Scenario")
        ax.set_ylabel("Predicted EC state")
        ax.set_xticks([0, 1])
        ax.set_xticklabels(cur.sort_values("x")["comparison_scenario"].drop_duplicates().tolist())
        ax.set_yticks(range(len(STATE_ORDER)))
        ax.set_yticklabels(STATE_ORDER)
        ax.set_xlim(-0.55, 1.55)
        ax.set_ylim(-0.55, len(STATE_ORDER) - 0.45)
        ax.grid(True, color="#e6e6e6", linewidth=1)
        ax.set_axisbelow(True)
    cbar = fig.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap=cmap),
        ax=axes,
        shrink=0.82,
        pad=0.03,
    )
    cbar.set_label("Frequency among all ECs")
    fig.savefig(OUT / "ec_state_frequency_among_all_ecs_pairwise_dotplot.pdf", bbox_inches="tight")
    fig.savefig(OUT / "ec_state_frequency_among_all_ecs_pairwise_dotplot.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def write_html():
    html_path = SRC / "mintflow_met_hi_fib_perturbation_ec_report.html"
    html = html_path.read_text()
    insertion = (
        "<h3>ec_state_frequency_among_all_ecs_pairwise_dotplot.png</h3>"
        "<img src='ec_state_frequency_among_all_ecs_pairwise_dotplot.png' style='max-width: 100%;'>"
    )
    marker = (
        "<h3>ec_state_frequency_among_all_ecs_pairwise.png</h3>"
        "<img src='ec_state_frequency_among_all_ecs_pairwise.png' style='max-width: 100%;'>"
    )
    if marker in html:
        html = html.replace(marker, marker + insertion, 1)
    else:
        html = re.sub(r"(<h2>Plots</h2>.*?</p>)", r"\1" + insertion, html, count=1, flags=re.S)
    highres_links = (
        "<h3>High-resolution tissue map PDFs</h3>"
        "<ul>"
        "<li>tissue_maps_celltype_perturbations_highres.pdf</li>"
        "<li>tissue_maps_ec_met_hi_probability_highres.pdf</li>"
        "<li>tissue_maps_met_hi_ec_met_hi_probability_highres.pdf</li>"
        "<li>tissue_maps_met_int_ec_met_hi_probability_highres.pdf</li>"
        "</ul>"
    )
    html = html.replace("<h2>Run Parameters</h2>", highres_links + "<h2>Run Parameters</h2>", 1)
    html = html.replace(
        "<title>MintFlow Met_hi_Fib Perturbation EC Analysis</title>",
        "<title>MintFlow Met_hi_Fib Perturbation EC Analysis with Frequency Dot Plots</title>",
    )
    html = html.replace(
        "<h1>MintFlow Met_hi_Fib in silico perturbation: EC response</h1>",
        "<h1>MintFlow Met_hi_Fib in silico perturbation: EC response with frequency dot plots</h1>",
    )
    (OUT / "mintflow_met_hi_fib_perturbation_ec_report_with_frequency_dotplots.html").write_text(html)


def main():
    copy_assets()
    make_pairwise_dotplot()
    write_html()
    print(OUT / "mintflow_met_hi_fib_perturbation_ec_report_with_frequency_dotplots.html")


if __name__ == "__main__":
    main()
