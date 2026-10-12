# Reproduction of recovered workflows (2026-10-12)

These additions recover numerical stages that were missing from the public repository. They are portable implementations of the original calculations, tested against separately saved project outputs. Private input objects and manual gate exports are required; the figure source-data tables do not replace these inputs.

| Added workflow | Reference comparison | Result |
|---|---|---|
| First-cohort area QC and expression transformation | 46,305 retained cells, 48 features in the saved Steinbock object | Cofactor-1 arcsinh, min–max and global z-score assays agree; largest absolute error 1.31e-14 |
| CD68/CD3 gate merge | Before/after first-cohort objects and 32 manual gate exports | All 45,340 cell labels match; all six assays unchanged |
| Histone gate merge | Before/after post-Xenium IMC objects and seven manual gate exports | All 6,218 labels match; 1,312 unlabeled Histone-positive cells become Stromal, 948 remain unlabeled; all five assays unchanged |
| Purity and silhouette scoring, k selection | December 2023 workspace, exported workbook and final object | All 32 rows and 320 rounded score values match exactly; all 45,340 selected cluster labels match |
| CX5 summary and statistics | Normalized workbook and archived rendered R Markdown report | All 12 summary rows and 10 tests match at the report's three-significant-digit display precision |

The preprocessing comparison starts with the retained-cell Steinbock object. It verifies the expression calculations and retained area range (20–250 inclusive), but does not independently reconstruct raw image segmentation or count excluded cells from raw Steinbock files. The Histone result belongs to the post-Xenium IMC cohort; it is distinct from the first cohort.

## Run the recovered stages

From the repository root, with the required R packages installed:

```sh
Rscript data_preprocessing/imc/recovered_workflows/preprocess_first_cohort.R unfiltered_spe.rds prepared_spe.rds
Rscript data_preprocessing/imc/recovered_workflows/merge_cd68_cd3_gates.R before.rds cd68_cd3_gate_exports after.rds
Rscript data_preprocessing/imc/recovered_workflows/merge_histone_gates.R before.rds histone_gate_exports after.rds
Rscript data_preprocessing/imc/recovered_workflows/score_purity_k_grid.R clustered_spe.rds data_preprocessing/imc/recovered_workflows/purity_markers_2023.txt purity
Rscript data_preprocessing/functional/recovered_workflows/analyze_cx5_normalized.R normalized.xlsx output
```

The merge scripts require RDS gate exports. The immune merge reads each export's cytomapper_CellLabel and limits reassignment to CD45-positive cells. The Histone merge uses cell IDs in exports containing only Histone-positive cells and limits reassignment to unlabeled cells. Manual gate placement itself is an upstream input, not an automated or revalidated decision. Source inputs are not modified.

The original immune script used an unfiltered file listing that would also try to read its own TXT file as RDS. The recovered script limits the listing to RDS files. Set-based cell matching also avoids redundant concatenation and preserves the verified labels. Explicit guards reject missing/unknown IDs and ambiguous inputs. No patient mappings or workstation paths are embedded.

## Purity interpretation

The recovered history uses k values 5, 20, 50, 100, 150, 200, 250 and 500 within four coarse compartments. Scoring uses 21 markers in the cell-type-standardized z_c_rescaled assay. neighborPurity uses k=50 and its historical weighted default; approxSilhouette uses centroid-based approximation. Selection maximizes mean purity across cells plus mean approximate silhouette across cells after rounding scores to three decimals. The source centers both terms without scaling; centering does not change their summed maximum.

The original exported headings mislabel columns 5/6 and 10/11. Those columns actually contain the mean and median of per-cluster means. The new scores CSV uses accurate names; the legacy CSV intentionally preserves the archived order, headings and rounding for comparison.

Despite the filename ending in 100_runs, the recovered clustering loop executes one seed (1902). Its clustering input is 11 pathway summaries, whereas its evaluation input is 21 marker values. The scorer accepts the archived Rphenograph_cluster_<k> columns and does not retrain Rphenograph. This archive selects k=500 for all four compartments, including unsplit CD45-positive cells. It does not independently establish the manuscript's later macrophage k=100 selection. Repeated-seed clustering and that final subset-specific choice remain unverified.

## CX5 interpretation

The recovered function starts from the Normalized_Data sheet with Channel, Condition and numeric Exp columns. The source analyzes donor/experiment-level percentages against 100 using two-sided one-sample t-tests and BH correction separately within each channel. Each condition has four observations except 2DG, which has three. All ten BH-adjusted comparisons are nonsignificant in this archived workbook. The original plot code separately calls unadjusted two-sample t-tests; its significance annotations are not the BH-adjusted table. These additions reproduce the numerical table only and do not relabel the original figure's statistical evidence.

Raw image feature extraction, DAPI normalization, technical-replicate aggregation and donor-to-experiment mapping have not been reconstructed from this workbook alone.

## Repeat the reference comparisons

Copy input_manifest.example.json to a private location and fill in the paths. The manifest requires the saved import object, before/after gating objects, gate folders, the December purity workspace/workbook/final object, and the normalized CX5 workbook. Run:

``sh
Rscript data_preprocessing/verification/reproduce_saved_outputs.R /private/input_manifest.json . /private/results.json
```

The script fails if any comparison fails and writes an aggregate results JSON plus numerical output tables beside it. Keep these outputs private if derived data are not cleared for sharing. The checked-in verification results contain aggregate evidence only. Source fingerprints are listed in RECOVERED_SOURCE_INDEX.md; runtime package versions are recorded separately in the validation metadata. The project's existing renv.lock is preserved and was not claimed to have been fully restored.

Other audit candidates—interactive gating apps, malformed alpha-SMA merge template, original six-donor Xenium annotation/ResolVI inputs, image registration transfer and Galaxy preprocessing—are not certified by these comparisons and have not been added as reproduced workflows.

## Windows runtime limitation observed during validation

All nine numerical/reference comparisons completed and passed on R 4.6.1 with the package versions recorded in runtime_versions.json. The R process then returned Windows access-violation code 0xC0000005. Loading igraph alone reproduces the shutdown fault; it also occurs with R 4.4.1/4.4.3 and the archived igraph 2.1.1 build. Base R and packages not importing igraph exit normally. The preserved numerical result files were written before shutdown, but a clean exit of the combined verification process is not certified. This environment fault has not been resolved by installing newer R and is explicitly recorded in verified_results.json.

Standalone CLI checks for preprocessing (R 4.6.1) and Histone merging (R 4.4.3 with Bioconductor 3.19) returned exit code zero and wrote their output objects. Old serialized objects can additionally require the legacy GenomeInfoDb namespace when saving with newer Bioconductor; use the original study environment or install the matching legacy namespace.

Histone merging and saving also passed on R 4.6.1 after installing GenomeInfoDb 1.40.1, GenomeInfoDbData 1.2.12 and UCSC.utils 1.0.0 from Bioconductor 3.19 for legacy serialized-object support. These R-only packages supplement the modern environment; the original renv.lock remains the recommended specification for historical sources.
