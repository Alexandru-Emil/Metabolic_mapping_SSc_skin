# Recovered workflow source index

Added on 2026-10-12 after numerical comparisons with preserved project outputs. These files are parameterized implementations of the recovered numerical stages, rather than byte-for-byte copies of the historical scripts. Original source SHA-256 fingerprints refer to untouched local originals. The import source credits Yi-Nan Li and the Bodenmiller Group ImagingWorkshop2023 workflow; the CX5 report is authored by Veda.

| Repository workflow | Original source relative to project folder | Original SHA-256 |
|---|---|---|
| [imc/recovered_workflows/preprocess_first_cohort.R](imc/recovered_workflows/preprocess_first_cohort.R) | K project: SOP and codes/20230516_1-read-steinbock.Rmd | `e69152a286c85c3b69c34ea84bb14e1401eeb8f260834670302a75e627618c02` |
| [imc/recovered_workflows/merge_cd68_cd3_gates.R](imc/recovered_workflows/merge_cd68_cd3_gates.R) | Sciebo project: CD68_CD3_gates/All gates/Script to merge and subset CD45 gate.txt | `01aa5f91e93c1533261e48ae71117925b6a25b800f6c4e4642f1acc2fa9f2f6f` |
| [imc/recovered_workflows/merge_histone_gates.R](imc/recovered_workflows/merge_histone_gates.R) | K project: SOP and codes/histone merged.R | `d995e86e74a926f6569a17a362153d417af1a7733c8ecd8175631b3cbae2089b` |
| [imc/recovered_workflows/score_purity_k_grid.R](imc/recovered_workflows/score_purity_k_grid.R) | K project: Metabolic panel/Results/Final figures/Purity run 5.12.2023/.Rhistory | `c07b7be31e4df9f912c5c681c21e85d5936df27aa70b90998d16d40961073666` |
| [functional/recovered_workflows/analyze_cx5_normalized.R](functional/recovered_workflows/analyze_cx5_normalized.R) | Sciebo project: Revision/RMD/SSc_fibroblast_analysis_NEW.Rmd | `4278d9cee82c14606413d8861d385d0fff315fd44274fa5af1b48fbc7c103f7c` |

The purity marker names were recovered from the archived workspace in the original z_total order. The full histories, source patient mapping, workstation paths, manual gate exports and input objects are not published. See [verification/README.md](verification/README.md) for comparisons, runtime metadata and limitations.
