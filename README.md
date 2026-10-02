# Metabolic mapping of SSc skin

Analysis code for **Single-cell mapping of the metabolic landscape of skin fibrosis in systemic sclerosis**, by Devakumar, Li, Filla, Rius Rigau, Györfi, Zhi, Zhu, Qi, Tümerdem, Neelagar, Wang, Liang, Bergmann, Schett, Distler and Matei.

The figure reference is the merged manuscript and figures in the 11 August 2026 revision. Numerical source data are distributed as a separate source-data archive.

## Organization

- [`data_preprocessing`](data_preprocessing/README.md): original IMC, Xenium, shared-modality, bulk RNA-seq, functional and MintFlow workflows, organized by purpose. Preparation notebooks and methods are retained; changes concern paths, source grouping, public donor aliases and dependency handling. [Source index](data_preprocessing/SOURCE_INDEX.md).
- [`figure_regeneration`](figure_regeneration/README.md): one script for each figure panel, named to match its source-data table. Numerical rendering starts with deposited CSV tables. [Panel index](figure_regeneration/PANEL_INDEX.md).

## Regenerate numerical panels

Use R 4.4.1 and the package versions in [`figure_regeneration/environment.lock`](figure_regeneration/environment.lock). Install dependencies into an environment you manage before running; the scripts do not install packages. Statistical summaries and tests are calculated from the supplied observations where needed to reproduce the displays.

Unzip `Source_data_repository_20261002.zip` into an external directory. Its SHA-256 is:

```text
0ca3798b4e8454c40b34445a9261ac2af4c327beb4909e570f217bf62f874640
```

From the repository root, run:

```sh
Rscript --vanilla figure_regeneration/verify_source_data.R /path/to/source_data
Rscript --vanilla figure_regeneration/render_all.R /path/to/source_data /path/to/rendered_panels
```

For one panel:

```sh
Rscript --vanilla figure_regeneration/Figure_5/Panel_E.R /path/to/source_data /path/to/rendered_panels
```

Windows paths containing spaces must be quoted. Outputs are PDF and PNG files arranged by figure, plus a table-verification report, rendered-panel index and R session information. Keep data and generated outputs outside the code checkout.

## Reproduction scope

The table-based scripts cover 146 numerical panels (567,421 source-data rows). Validation checks require the exact source-table checksums, column schemas and row counts, including the final six-donor validation cohort, 622 metabolically high fibroblasts and the eight measurements per collagen condition. Recorded F/M categories accompany individual IMC and shared IMC–Xenium records where available; categories are not inferred for other experiments.

Original plotting code and components are used with table-loading adapters. The scripts regenerate numerical plots; final multi-panel manuscript layouts, edited legends, clinical annotation strips not represented in the deposited tables, schematics and representative microscopy images require their original assets. A complete, pixel-identical recreation of every assembled manuscript page is therefore outside the numerical-table run.

S9D uses the reference-projection coordinates regenerated from the original PCA/UMAP/anchor/MapQuery workflow on 2 October 2026. These are newly computed coordinates rather than recovered original saved coordinates. The table-rendering script draws those deposited coordinates without recomputing an embedding.

Figure 5G's supplied table contains raw CoreMatrisome AUCell means; the 11 August figure's negative colour legend requires a different or standardized display scale whose provenance has not been established. The available raw-mean plot runs, but its exact agreement with that panel is not certified.

The preparation workflows were checked for source syntax and organization. The raw-data preparation pipeline and model training were not rerun as part of this code-release validation. [Verification details](figure_regeneration/VERIFICATION.md).

## Provenance and attribution

The source review includes [IMC_metabolic_AEM, commit 3e1a84b](https://github.com/Alexandru-Emil/IMC_metabolic_AEM/tree/3e1a84b7f0be9340604765f1a02255435b4aabcd) and project-local analysis revisions. Original and public file fingerprints are recorded in the source index. Earlier seven-donor workflows are retained in explicitly named comparative folders; they are separate from the final six-donor manuscript figure inputs.

The original IMC import workflow acknowledges Nils Eling and the [Bodenmiller Group Imaging Workshop 2023](https://github.com/BodenmillerGroup/ImagingWorkshop2023). Package citations and third-party notices remain applicable. This repository does not contain patient-identifying source files, original clinical workbooks, manuscript documents or generated reports.
