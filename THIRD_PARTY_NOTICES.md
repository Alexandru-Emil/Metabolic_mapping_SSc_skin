# Third-party notices

## Bodenmiller Group Imaging Workshop 2023

The IMC import notebook [`data_preprocessing/imc/workflows/20230412_1-read-steinbock.Rmd`](data_preprocessing/imc/workflows/20230412_1-read-steinbock.Rmd) contains text and code adapted by Yi-Nan Li from Nils Eling's material for the [Multiplexed Tissue Imaging Workshop 2023](https://github.com/BodenmillerGroup/ImagingWorkshop2023). Its original acknowledgement is retained.

Upstream copyright: **Copyright (c) 2022 BodenmillerGroup**. Upstream licence: **MIT**. The complete, unchanged notice is retained in [`ImagingWorkshop2023_MIT_LICENSE.txt`](data_preprocessing/third_party_licenses/ImagingWorkshop2023_MIT_LICENSE.txt).

The licence was verified at upstream commit `970381c43efc8f6cfd1ff6014950613edbdd1f50` on 2 October 2026. This identifies the licence checked, rather than asserting that the 2023 adaptation used that later code revision.

## Bodenmiller Group IMCDataAnalysis

Project IMC workflows also refer to the [IMCDataAnalysis workflow](https://bodenmillergroup.github.io/IMCDataAnalysis/), including its image- and cell-level quality-control guidance. Its source repository is [BodenmillerGroup/IMCDataAnalysis](https://github.com/BodenmillerGroup/IMCDataAnalysis). The complete licence is retained to accompany any adapted workflow code and documentation.

Upstream copyright: **Copyright (c) 2020 BodenmillerGroup**. Upstream licence: **MIT**. The complete, unchanged notice is retained in [`IMCDataAnalysis_MIT_LICENSE.txt`](data_preprocessing/third_party_licenses/IMCDataAnalysis_MIT_LICENSE.txt).

The licence was verified at upstream commit `89bfb3f0f2d5732334bb26a8aa32847bbc67ff81` on 2 October 2026.

The workflow's recommended scientific citation is Windhager, J., Zanotelli, V. R. T., Schulz, D. et al. *An end-to-end workflow for multiplexed image processing and analysis*. Nature Protocols (2023). [https://doi.org/10.1038/s41596-023-00881-0](https://doi.org/10.1038/s41596-023-00881-0).

## Licence scope and dependencies

The root [`LICENSE`](LICENSE) grants the MIT licence for this repository's project code and accompanying documentation and retains the above upstream copyright notices. It does not replace the licences of external R, Python or other dependencies. Those packages remain subject to their own licences and scientific citation requirements; dependency versions are recorded in the environment files.

The MIT code licence does not grant rights to separately deposited clinical or experimental data, microscopy images, manuscript assets or third-party datasets. Their repository records specify the applicable access conditions and reuse terms.

## Public sample labels

The historical import notebook uses `ExampleHealthy1`, `ExampleHealthy2`, `ExampleSSc1` and `ExampleSSc2` in place of its original acquisition and patient labels. These are illustrative historical labels and do not establish a mapping to the final deposited cohorts. The original source fingerprint remains in the preparation-source index.

Legacy healthy-line labels were removed from Seahorse explanatory text. EC source codes and the historical additional Xenium sample label are supplied through private runtime configuration when needed; numerical analysis and figure-table inputs are unchanged. Earlier Git commits have not been rewritten.
