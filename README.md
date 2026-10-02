# Metabolic mapping of SSc skin

Analysis code for **Single-cell mapping of the metabolic landscape of skin fibrosis in systemic sclerosis**, by Devakumar, Li, Filla, Rius Rigau, Györfi, Zhi, Zhu, Qi, Tümerdem, Neelagar, Wang, Liang, Bergmann, Schett, Distler and Matei.

Archived code version 1.0.0: [10.5281/zenodo.23107219](https://doi.org/10.5281/zenodo.23107219). Figure source data: [10.5281/zenodo.23107366](https://doi.org/10.5281/zenodo.23107366).

## Organization

- [`data_preprocessing`](data_preprocessing/README.md): original IMC, Xenium, shared-modality, bulk RNA-seq, functional and MintFlow workflows, organized by purpose. [Source index](data_preprocessing/SOURCE_INDEX.md).
- [`figure_regeneration`](figure_regeneration/README.md): one script for each figure panel, named to match its source-data table. [Panel index](figure_regeneration/PANEL_INDEX.md).

## Regenerate numerical panels

Use R 4.4.1 and the package versions in [`figure_regeneration/environment.lock`](figure_regeneration/environment.lock). 

Unzip `Source_data_repository_20261002.zip` into an external directory. Its SHA-256 is:

```text
0b26712b2ce350852e102f751cc2a8206c41bd7cc415aa4b21467a59349ba7e6
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

## Provenance and attribution

The original IMC import workflow acknowledges Nils Eling and the [Bodenmiller Group Imaging Workshop 2023](https://github.com/BodenmillerGroup/ImagingWorkshop2023). Package citations and third-party notices remain applicable. Clinical source workbooks, manuscript documents and generated reports are distributed outside this code repository.

## Licence and citation

This repository's code and accompanying documentation are available under the [MIT licence](LICENSE). Copyright notices and complete licences for adapted Bodenmiller Group material are preserved in [third-party notices](THIRD_PARTY_NOTICES.md). External dependencies and separately deposited data retain their own licences and access conditions.
