# Reproduction verification

The final table archive is `Source_data_repository_20261002.zip` (SHA-256 `0ca3798b4e8454c40b34445a9261ac2af4c327beb4909e570f217bf62f874640`). The public renderer is tested against its unpacked contents, without raw objects, cached model files, Excel or embedded Prism files as plotting inputs.

## Checks

- All 146 numerical panel CSV files match the recorded SHA-256 checksums, expected column order and row counts: 567,421 rows in total.
- Each named numerical panel script runs and exports a nonempty PDF and PNG. The full run writes `rendered_panels.csv`, `verified_source_tables.csv` and `sessionInfo.txt` outside the code repository.
- The renderer builds each panel layout once with fixed device settings and reuses its grobs for PDF and PNG. The random seed is reset before construction and drawing, including jittered observations and repelled labels.
- The full run executes the same named scripts as the single-panel commands, each in a separate R session, to prevent graphics and font state carrying over between panels.
- Coordinate columns contain finite values. Frequencies, percentages and niche fractions agree with their stored counts and denominators where these are supplied.
- F/M records are complete in the individual IMC/shared-modality tables that include recorded categories. Other experiments do not receive inferred categories.
- The validation fibroblast UMAP includes 1,155 cells from six donors, with 622 Met_hi_Fib cells. S9D includes 20,657 reference and 1,155 projected validation cells. The collagen comparison includes eight measurements in each of three conditions.
- Preparation R/Rmd source code and Python files pass syntax checks. No preprocessing, annotation or model-training rerun is implied by these source checks.

## Figure-reference comparison

The source-table checks establish which numerical inputs are used and whether the named scripts execute. They are distinct from a pixel comparison with an assembled manuscript figure.

Original six-donor IMC plot expressions are retained for 5B/C/E, S9A/B and S14A/B, including their colour limits, paired statistics and export dimensions. Original functions are retained for signal/noise, GSEA and metabolic/ECM score displays. Other panels use the original plotting components with table bindings. Prism-derived functional plots are recreated from numerical observations using the original plot type and experimental design.

Figure 5D uses a taller export to keep all gene labels separate, with the original marker/colour scales and the manuscript's horizontal group labels and legend order. Figure 5F retains the original Sankey calculations, with the manuscript's node palette, pale flows and external labels. Numerical inputs are unchanged by these presentation settings.

Figure 5G's recovered table contains nonnegative CoreMatrisome AUCell means. The 11 August figure displays negative values in its ECM colour legend. The supplied adapter draws the recovered raw means; the original figure's displayed colour transformation or alternative score has not been established from the located sources. An exact Figure 5G match is therefore not certified.

S9D draws newly regenerated reference-projection coordinates from the original PCA/UMAP/anchor/MapQuery method. These coordinates replace a missing saved projection and are not claimed to be the original unsaved embedding.

Final page layouts, manually edited labels/legends, clinical annotation strips absent from the panel tables, schematics and microscopy images are supplied by the original manuscript assets. The 19 image/schematic wrapper scripts require original PNG assets and are separate from the numerical-table run.
