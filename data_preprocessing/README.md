# Data preparation workflows

The source files preserve the project analyses and their alternative versions. They are organized by modality and purpose, with the original path and file fingerprint recorded in [SOURCE_INDEX.md](SOURCE_INDEX.md). R Markdown narrative, chunks and evaluation flags are retained.

| Folder | Purpose |
|---|---|
| `imc` | Steinbock import, rescaling, cell/marker QC, clustering, scores, frequencies, distances, neighborhoods, interactions and clinical associations |
| `xenium` | Spatial object and sample mapping, cell filtering, annotation and ResolVI workflows |
| `shared` | Matched IMC–Xenium labels, metabolic and ECM scores, cohort validation and reference projection |
| `rnaseq` | Bulk RNA-seq QC, batch correction, edgeR/DESeq2 comparisons, selected genes and enrichment |
| `functional` | Lactate, Seahorse and crossed-design experimental analysis |
| `mintflow` | Model input export, training, perturbation and cached-output analysis |
| `common` | Shared input configuration and original plotting helpers |

`workflows` contains the project-local sources. `additional_source_workflows` retains distinct scientific sources from the existing GitHub repository. `comparative_cohort_workflows` explicitly preserves seven-donor reanalyses for traceability. Hash prefixes distinguish source versions with the same original filename.

The main validation figure code is `shared/workflows/six_donor_imc_only_figure_panels.Rmd`. Its cohort is six donors. The figure-rendering layer uses the final deposited tables and never selects a seven-donor variant.

Set `METABOLIC_INPUT_DIR` to an external directory containing the project-relative input hierarchy; set `METABOLIC_OUTPUT_DIR` for analysis outputs. Run from the repository root or set the R Markdown knitting root accordingly. Inputs such as Steinbock files, SpatialExperiment/Seurat objects, raw counts, matrisome definitions and model checkpoints are required by their individual workflows and are not substituted by the small figure tables.

The legacy source workflows retain their original internal ordering and object dependencies. Choose the intended source version and required upstream objects using the source index. They are not an automatically chained pipeline. Existing methods, thresholds, transformations, donor exclusions and alternative analyses have not been redefined for this release.

The supplied `renv.lock` records the original IMC project's dependency specification. The separately tested table-rendering environment is recorded under `figure_regeneration`.

## Legacy source-code configuration

The historical Seahorse preparation workflow reads three EC source codes from `METABOLIC_EC_SOURCE_CODES`, as a comma-separated string in the original EC-line 1, 2, 3 order. The default is `EC1,EC2,EC3` for an aliased input workbook. Use a private environment value when working with the original workbook; source donor codes are not embedded in the public script. This changes input configuration only, preserving the original factor order and statistical design.

The historical ResolVI input notebook reads its additional sample label from `METABOLIC_ADDITIONAL_XENIUM_SAMPLE`, defaulting to `AdditionalValidation1`. This additional sample is separate from the final six-donor deposited figure cohort. The excluded-validation QC workflows use descriptive variable names.
