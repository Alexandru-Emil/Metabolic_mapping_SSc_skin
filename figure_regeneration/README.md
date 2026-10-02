# Figure regeneration from deposited tables

Each `Figure_<id>/Panel_<letter>.R` reads `Figure_<id>/Panel_<letter>_values.csv`. These scripts use the shared renderer, original plotting expressions and panel-specific adapters. [PANEL_INDEX.md](PANEL_INDEX.md) links each plot to its source table and original workflow.

Run from the repository root, giving the external data and output directories as arguments:

```sh
Rscript --vanilla figure_regeneration/render_all.R DATA_DIR OUTPUT_DIR
Rscript --vanilla figure_regeneration/Figure_3/Panel_C.R DATA_DIR OUTPUT_DIR
```

The source-data directory must contain `File_manifest.csv`, `Cell_state_group_labels.csv` and all 146 panel CSV files. `verify_source_data.R` rejects missing, altered or structurally incompatible tables. The renderer uses stored plotting values and coordinates rather than recalculating clustering, scores, expression normalization, differential expression or model training.

Statistics and plot summaries needed to draw the original displays (e.g. paired tests, boxplot quartiles, mean/SEM and density estimates) are computed from the supplied numerical observations. The original statistical design is retained where applicable. The Seahorse endpoints use the original crossed design with EC line as a fixed effect and fibroblast donor as a random intercept.

The full run launches each named panel script in a separate R session to keep graphics and font settings independent between panels. Four panels run concurrently by default; set `METABOLIC_PLOT_WORKERS=1` for sequential execution on a computer with limited memory. Panel logs and session information accompany the exports.

Niche fractions and frequencies defined as missing because of insufficient cells retain missing values. They are not replaced with zero. Tables may include additional scores; adapters select the scores shown in the corresponding manuscript panel.

The 19 image/schematic panel scripts use original PNG assets instead of numerical tables:

```sh
Rscript --vanilla figure_regeneration/Figure_1/Panel_A.R ASSET_DIR OUTPUT_DIR
```

Supply `ASSET_DIR/Figure_1/Panel_A_image.png` (and equivalent names for other image panels). Mixed microscopy/quantification panels need the original microscopy images for final assembly even where their numerical plot is regenerated here.

`environment.lock` records package versions and their installation provenance from the tested R environment. No automatic dependency installation occurs. See [VERIFICATION.md](VERIFICATION.md) for the scope of verification and the S9D coordinate provenance.
