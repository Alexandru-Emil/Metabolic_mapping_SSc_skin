# Metabolic mapping of SSc skin

Analysis code for **Single-cell mapping of the metabolic landscape of skin fibrosis in systemic sclerosis**, by Devakumar, Li, Filla, Rius Rigau, Györfi, Zhi, Zhu, Qi, Tümerdem, Neelagar, Wang, Liang, Bergmann, Schett, Distler and Matei.

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


## Provenance and attribution

The original IMC import workflow acknowledges Nils Eling and the [Bodenmiller Group Imaging Workshop 2023](https://github.com/BodenmillerGroup/ImagingWorkshop2023). Package citations and third-party notices remain applicable. This repository does not contain patient-identifying source files, original clinical workbooks, manuscript documents or generated reports.
