# Panel-to-table and plotting-source index

Run all numerical panels with `Rscript figure_regeneration/render_all.R DATA_DIR OUTPUT_DIR` from the repository root. Each named script reads its matching deposited table and the shared label file. Panels 6D/E also read each other's table to coordinate their axes and legends. The registry explicitly fixes the column schema, row count and plot specification.

Original plotting functions and expressions are retained in `common/original_plot_code.R`; adapters replace analysis-object bindings with deposited tables. Other adapters retain the original plotting components (e.g. dittoDimPlot, ggradar, ggsankey, distance-density and interaction plots). Embedded Prism plots are recreated from their exported numerical values. Manuscript page assembly and image pixels require the original assets.

| Figure panel | Script | Source table | Rows | Original plotting source/profile |
|---|---|---|---:|---|
| 10A | [Figure_10/Panel_A.R](Figure_10/Panel_A.R) | `Figure_10/Panel_A_values.csv` | 10 | imc/workflows/first_cohort_mrss_common.R; imc/workflows/Regression_clinicaldata/Clinical_correlations_dotplots.Rmd |
| 10B | [Figure_10/Panel_B.R](Figure_10/Panel_B.R) | `Figure_10/Panel_B_values.csv` | 9 | imc/workflows/first_cohort_mrss_common.R; imc/workflows/Regression_clinicaldata/Clinical_correlations_dotplots.Rmd |
| 10C | [Figure_10/Panel_C.R](Figure_10/Panel_C.R) | `Figure_10/Panel_C_values.csv` | 9 | imc/workflows/first_cohort_mrss_common.R; imc/workflows/Regression_clinicaldata/Clinical_correlations_dotplots.Rmd |
| 10D | [Figure_10/Panel_D.R](Figure_10/Panel_D.R) | `Figure_10/Panel_D_values.csv` | 9 | imc/workflows/first_cohort_mrss_common.R; imc/workflows/Regression_clinicaldata/Clinical_correlations_dotplots.Rmd |
| 10E | [Figure_10/Panel_E.R](Figure_10/Panel_E.R) | `Figure_10/Panel_E_values.csv` | 9 | imc/workflows/first_cohort_mrss_common.R; imc/workflows/Regression_clinicaldata/Clinical_correlations_dotplots.Rmd |
| 10F | [Figure_10/Panel_F.R](Figure_10/Panel_F.R) | `Figure_10/Panel_F_values.csv` | 9 | imc/workflows/first_cohort_mrss_common.R; imc/workflows/Regression_clinicaldata/Clinical_correlations_dotplots.Rmd |
| 10G | [Figure_10/Panel_G.R](Figure_10/Panel_G.R) | `Figure_10/Panel_G_values.csv` | 9 | imc/workflows/first_cohort_mrss_common.R; imc/workflows/Regression_clinicaldata/Clinical_correlations_dotplots.Rmd |
| 10H | [Figure_10/Panel_H.R](Figure_10/Panel_H.R) | `Figure_10/Panel_H_values.csv` | 9 | imc/workflows/first_cohort_mrss_common.R; imc/workflows/Regression_clinicaldata/Clinical_correlations_dotplots.Rmd |
| 2A | [Figure_2/Panel_A.R](Figure_2/Panel_A.R) | `Figure_2/Panel_A_values.csv` | 336 | imc/workflows/Met_expr/Met_expr.Rmd; imc/workflows/20231212_FibClustering.Rmd |
| 2B | [Figure_2/Panel_B.R](Figure_2/Panel_B.R) | `Figure_2/Panel_B_values.csv` | 16 | imc/workflows/Met_expr/Met_expr.Rmd |
| 2C | [Figure_2/Panel_C.R](Figure_2/Panel_C.R) | `Figure_2/Panel_C_values.csv` | 64 | imc/workflows/Freq/Freq_clust.Rmd; imc/workflows/Met_expr/Met_expr.Rmd |
| 3A | [Figure_3/Panel_A.R](Figure_3/Panel_A.R) | `Figure_3/Panel_A_values.csv` | 20,657 | imc/workflows/20231212_FibClustering.Rmd; shared/workflows/IMC_heatmap_sankey_new.Rmd |
| 3B | [Figure_3/Panel_B.R](Figure_3/Panel_B.R) | `Figure_3/Panel_B_values.csv` | 20,657 | imc/workflows/20231212_FibClustering.Rmd |
| 3C | [Figure_3/Panel_C.R](Figure_3/Panel_C.R) | `Figure_3/Panel_C_values.csv` | 168 | imc/workflows/Met_expr/Met_expr.Rmd; imc/workflows/20231212_FibClustering.Rmd |
| 3D | [Figure_3/Panel_D.R](Figure_3/Panel_D.R) | `Figure_3/Panel_D_values.csv` | 24 | imc/workflows/Met_expr/Met_expr.Rmd; imc/workflows/20231212_FibClustering.Rmd |
| 4A | [Figure_4/Panel_A.R](Figure_4/Panel_A.R) | `Figure_4/Panel_A_values.csv` | 20,657 | imc/workflows/20231212_FibClustering.Rmd; shared/workflows/IMC_heatmap_sankey_new.Rmd |
| 4B | [Figure_4/Panel_B.R](Figure_4/Panel_B.R) | `Figure_4/Panel_B_values.csv` | 24 | imc/workflows/Freq/Freq_clust.Rmd |
| 4C | [Figure_4/Panel_C.R](Figure_4/Panel_C.R) | `Figure_4/Panel_C_values.csv` | 16 | imc/workflows/Freq/Freq_clust.Rmd; imc/workflows/Met_expr/Met_expr.Rmd |
| 4D | [Figure_4/Panel_D.R](Figure_4/Panel_D.R) | `Figure_4/Panel_D_values.csv` | 16 | imc/workflows/Freq/Freq_clust.Rmd; imc/workflows/Met_expr/Met_expr.Rmd |
| 5B | [Figure_5/Panel_B.R](Figure_5/Panel_B.R) | `Figure_5/Panel_B_values.csv` | 1,155 | shared/workflows/six_donor_imc_only_figure_panels.Rmd |
| 5C | [Figure_5/Panel_C.R](Figure_5/Panel_C.R) | `Figure_5/Panel_C_values.csv` | 1,155 | shared/workflows/six_donor_imc_only_figure_panels.Rmd |
| 5D | [Figure_5/Panel_D.R](Figure_5/Panel_D.R) | `Figure_5/Panel_D_values.csv` | 44 | shared/workflows/IMC_heatmap_sankey_new.Rmd |
| 5E | [Figure_5/Panel_E.R](Figure_5/Panel_E.R) | `Figure_5/Panel_E_values.csv` | 36 | shared/workflows/six_donor_imc_only_figure_panels.Rmd |
| 5F | [Figure_5/Panel_F.R](Figure_5/Panel_F.R) | `Figure_5/Panel_F_values.csv` | 16 | shared/workflows/IMC_heatmap_sankey_new.Rmd |
| 5G | [Figure_5/Panel_G.R](Figure_5/Panel_G.R) | `Figure_5/Panel_G_values.csv` | 64 | shared/workflows/Xenium_fibroblast_ECM_scores_by_IMC_metfiblabel_report.Rmd |
| 6B | [Figure_6/Panel_B.R](Figure_6/Panel_B.R) | `Figure_6/Panel_B_values.csv` | 2,140 | rnaseq/workflows/15_all_feature_strict_deg_figure_heatmaps_edgeR_exact.R |
| 6C | [Figure_6/Panel_C.R](Figure_6/Panel_C.R) | `Figure_6/Panel_C_values.csv` | 1,410 | rnaseq/workflows/15_all_feature_strict_deg_figure_heatmaps_edgeR_exact.R |
| 6D | [Figure_6/Panel_D.R](Figure_6/Panel_D.R) | `Figure_6/Panel_D_values.csv` | 11 | rnaseq/workflows/11_selected_2dg_ldha_gsea_dotplots_edgeR_exact_merged_metabolism.R |
| 6E | [Figure_6/Panel_E.R](Figure_6/Panel_E.R) | `Figure_6/Panel_E_values.csv` | 11 | rnaseq/workflows/11_selected_2dg_ldha_gsea_dotplots_edgeR_exact_merged_metabolism.R |
| 7A | [Figure_7/Panel_A.R](Figure_7/Panel_A.R) | `Figure_7/Panel_A_values.csv` | 16 | imc/workflows/Interaction/Expr_hasneigh_plots.Rmd |
| 7B | [Figure_7/Panel_B.R](Figure_7/Panel_B.R) | `Figure_7/Panel_B_values.csv` | 15 | imc/workflows/Interaction/Expr_hasneigh_plots.Rmd |
| 7C | [Figure_7/Panel_C.R](Figure_7/Panel_C.R) | `Figure_7/Panel_C_values.csv` | 4,373 | imc/workflows/Dist_ridge/Dist_to_cell_types_density_heatmap.Rmd |
| 7D | [Figure_7/Panel_D.R](Figure_7/Panel_D.R) | `Figure_7/Panel_D_values.csv` | 727 | imc/workflows/Dist_ridge/Dist_to_cell_types_density_heatmap.Rmd |
| 7E | [Figure_7/Panel_E.R](Figure_7/Panel_E.R) | `Figure_7/Panel_E_values.csv` | 20 | imc/workflows/Dist_ridge/Dist_to_cell_types_density_heatmap.Rmd |
| 7F | [Figure_7/Panel_F.R](Figure_7/Panel_F.R) | `Figure_7/Panel_F_values.csv` | 16 | imc/workflows/Dist_ridge/Dist_to_cell_types_density_heatmap.Rmd |
| 8A | [Figure_8/Panel_A.R](Figure_8/Panel_A.R) | `Figure_8/Panel_A_values.csv` | 190 | imc/workflows/Met_expr/Met_expr.Rmd; imc/workflows/20231212_FibClustering.Rmd |
| 8B | [Figure_8/Panel_B.R](Figure_8/Panel_B.R) | `Figure_8/Panel_B_values.csv` | 43 | imc/workflows/Interaction/Interaction_tests_plots_beautiful_plots.Rmd |
| 9A | [Figure_9/Panel_A.R](Figure_9/Panel_A.R) | `Figure_9/Panel_A_values.csv` | 13 | imc/workflows/Regression_clinicaldata/Differences_binary_clinical_outcome_dotplots.Rmd |
| 9B | [Figure_9/Panel_B.R](Figure_9/Panel_B.R) | `Figure_9/Panel_B_values.csv` | 9 | imc/workflows/Freq/Freq_clust.Rmd; imc/workflows/Met_expr/Met_expr.Rmd |
| 9C | [Figure_9/Panel_C.R](Figure_9/Panel_C.R) | `Figure_9/Panel_C_values.csv` | 9 | imc/workflows/Freq/Freq_clust.Rmd; imc/workflows/Met_expr/Met_expr.Rmd |
| 9D | [Figure_9/Panel_D.R](Figure_9/Panel_D.R) | `Figure_9/Panel_D_values.csv` | 9 | imc/workflows/Freq/Freq_clust.Rmd; imc/workflows/Met_expr/Met_expr.Rmd |
| 9E | [Figure_9/Panel_E.R](Figure_9/Panel_E.R) | `Figure_9/Panel_E_values.csv` | 9 | imc/workflows/Freq/Freq_clust.Rmd; imc/workflows/Met_expr/Met_expr.Rmd |
| 9F | [Figure_9/Panel_F.R](Figure_9/Panel_F.R) | `Figure_9/Panel_F_values.csv` | 9 | imc/workflows/Freq/Freq_clust.Rmd; imc/workflows/Met_expr/Met_expr.Rmd |
| 9G | [Figure_9/Panel_G.R](Figure_9/Panel_G.R) | `Figure_9/Panel_G_values.csv` | 9 | imc/workflows/Freq/Freq_clust.Rmd; imc/workflows/Met_expr/Met_expr.Rmd |
| 9H | [Figure_9/Panel_H.R](Figure_9/Panel_H.R) | `Figure_9/Panel_H_values.csv` | 9 | imc/workflows/Freq/Freq_clust.Rmd; imc/workflows/Met_expr/Met_expr.Rmd |
| S10A | [Figure_S10/Panel_A.R](Figure_S10/Panel_A.R) | `Figure_S10/Panel_A_values.csv` | 24 | Original embedded Prism/Excel numerical plots; bar/SEM adapter |
| S10B | [Figure_S10/Panel_B.R](Figure_S10/Panel_B.R) | `Figure_S10/Panel_B_values.csv` | 24 | Original embedded Prism/Excel numerical plots; bar/SEM adapter |
| S10C | [Figure_S10/Panel_C.R](Figure_S10/Panel_C.R) | `Figure_S10/Panel_C_values.csv` | 24 | Original embedded Prism/Excel numerical plots; bar/SEM adapter |
| S11A | [Figure_S11/Panel_A.R](Figure_S11/Panel_A.R) | `Figure_S11/Panel_A_values.csv` | 48 | Original embedded Prism/Excel numerical plots; bar/SEM adapter |
| S11B | [Figure_S11/Panel_B.R](Figure_S11/Panel_B.R) | `Figure_S11/Panel_B_values.csv` | 24 | Original embedded Prism/Excel numerical plots; bar/SEM adapter |
| S12unlettered | [Figure_S12/Panel_unlettered.R](Figure_S12/Panel_unlettered.R) | `Figure_S12/Panel_unlettered_values.csv` | 12 | Archived functional quantification; normalized bar/SEM adapter |
| S13A | [Figure_S13/Panel_A.R](Figure_S13/Panel_A.R) | `Figure_S13/Panel_A_values.csv` | 4,373 | imc/workflows/20231212_FibClustering.Rmd; shared/workflows/IMC_heatmap_sankey_new.Rmd |
| S13B | [Figure_S13/Panel_B.R](Figure_S13/Panel_B.R) | `Figure_S13/Panel_B_values.csv` | 4,373 | imc/workflows/20231212_FibClustering.Rmd |
| S13C | [Figure_S13/Panel_C.R](Figure_S13/Panel_C.R) | `Figure_S13/Panel_C_values.csv` | 105 | imc/workflows/Met_expr/Met_expr.Rmd; imc/workflows/20231212_FibClustering.Rmd |
| S13D | [Figure_S13/Panel_D.R](Figure_S13/Panel_D.R) | `Figure_S13/Panel_D_values.csv` | 4,373 | imc/workflows/Freq/Freq_clust.Rmd; imc/workflows/Met_expr/Met_expr.Rmd |
| S14A | [Figure_S14/Panel_A.R](Figure_S14/Panel_A.R) | `Figure_S14/Panel_A_values.csv` | 350 | shared/workflows/six_donor_imc_only_figure_panels.Rmd |
| S14B | [Figure_S14/Panel_B.R](Figure_S14/Panel_B.R) | `Figure_S14/Panel_B_values.csv` | 42 | shared/workflows/six_donor_imc_only_figure_panels.Rmd |
| S14C | [Figure_S14/Panel_C.R](Figure_S14/Panel_C.R) | `Figure_S14/Panel_C_values.csv` | 508 | shared/workflows/IMC_protein_and_scMetabolism_scores_combined.Rmd |
| S14D | [Figure_S14/Panel_D.R](Figure_S14/Panel_D.R) | `Figure_S14/Panel_D_values.csv` | 1,369 | imc/workflows/20231212_FibClustering.Rmd; shared/workflows/IMC_heatmap_sankey_new.Rmd |
| S14E | [Figure_S14/Panel_E.R](Figure_S14/Panel_E.R) | `Figure_S14/Panel_E_values.csv` | 72 | shared/workflows/IMC_heatmap_sankey_new.Rmd |
| S14F | [Figure_S14/Panel_F.R](Figure_S14/Panel_F.R) | `Figure_S14/Panel_F_values.csv` | 5,476 | shared/workflows/Xenium_EC_ECM_scores_report.Rmd |
| S14G | [Figure_S14/Panel_G.R](Figure_S14/Panel_G.R) | `Figure_S14/Panel_G_values.csv` | 2 | shared/workflows/IMC_heatmap_sankey_new.Rmd |
| S14H | [Figure_S14/Panel_H.R](Figure_S14/Panel_H.R) | `Figure_S14/Panel_H_values.csv` | 2,738 | shared/workflows/ACTA2_EC_scMetabolism_report.Rmd |
| S15B | [Figure_S15/Panel_B.R](Figure_S15/Panel_B.R) | `Figure_S15/Panel_B_values.csv` | 15 | imc/workflows/Freq/Freq_clust.Rmd |
| S15C | [Figure_S15/Panel_C.R](Figure_S15/Panel_C.R) | `Figure_S15/Panel_C_values.csv` | 16 | imc/workflows/Freq/Freq_clust.Rmd; imc/workflows/Met_expr/Met_expr.Rmd |
| S15D | [Figure_S15/Panel_D.R](Figure_S15/Panel_D.R) | `Figure_S15/Panel_D_values.csv` | 16 | imc/workflows/Freq/Freq_clust.Rmd; imc/workflows/Met_expr/Met_expr.Rmd |
| S16A | [Figure_S16/Panel_A.R](Figure_S16/Panel_A.R) | `Figure_S16/Panel_A_values.csv` | 727 | imc/workflows/20231212_FibClustering.Rmd; shared/workflows/IMC_heatmap_sankey_new.Rmd |
| S16B | [Figure_S16/Panel_B.R](Figure_S16/Panel_B.R) | `Figure_S16/Panel_B_values.csv` | 727 | imc/workflows/20231212_FibClustering.Rmd |
| S16C | [Figure_S16/Panel_C.R](Figure_S16/Panel_C.R) | `Figure_S16/Panel_C_values.csv` | 84 | imc/workflows/Met_expr/Met_expr.Rmd; imc/workflows/20231212_FibClustering.Rmd |
| S17B | [Figure_S17/Panel_B.R](Figure_S17/Panel_B.R) | `Figure_S17/Panel_B_values.csv` | 12 | imc/workflows/Freq/Freq_clust.Rmd |
| S17C | [Figure_S17/Panel_C.R](Figure_S17/Panel_C.R) | `Figure_S17/Panel_C_values.csv` | 16 | imc/workflows/Freq/Freq_clust.Rmd; imc/workflows/Met_expr/Met_expr.Rmd |
| S17D | [Figure_S17/Panel_D.R](Figure_S17/Panel_D.R) | `Figure_S17/Panel_D_values.csv` | 16 | imc/workflows/Freq/Freq_clust.Rmd; imc/workflows/Met_expr/Met_expr.Rmd |
| S18A | [Figure_S18/Panel_A.R](Figure_S18/Panel_A.R) | `Figure_S18/Panel_A_values.csv` | 45,340 | imc/workflows/20231212_FibClustering.Rmd; shared/workflows/IMC_heatmap_sankey_new.Rmd |
| S18B | [Figure_S18/Panel_B.R](Figure_S18/Panel_B.R) | `Figure_S18/Panel_B_values.csv` | 45,340 | imc/workflows/20231212_FibClustering.Rmd; shared/workflows/IMC_heatmap_sankey_new.Rmd |
| S19A | [Figure_S19/Panel_A.R](Figure_S19/Panel_A.R) | `Figure_S19/Panel_A_values.csv` | 20,657 | imc/workflows/Dist_ridge/Dist_to_cell_types_density_heatmap.Rmd |
| S19B | [Figure_S19/Panel_B.R](Figure_S19/Panel_B.R) | `Figure_S19/Panel_B_values.csv` | 20,657 | imc/workflows/Dist_ridge/Dist_to_cell_types_density_heatmap.Rmd |
| S19C | [Figure_S19/Panel_C.R](Figure_S19/Panel_C.R) | `Figure_S19/Panel_C_values.csv` | 20,657 | imc/workflows/Dist_ridge/Dist_to_cell_types_density_heatmap.Rmd |
| S2A | [Figure_S2/Panel_A.R](Figure_S2/Panel_A.R) | `Figure_S2/Panel_A_values.csv` | 73 | imc/workflows/compute_imc_signal_to_noise_qc.R |
| S2B | [Figure_S2/Panel_B.R](Figure_S2/Panel_B.R) | `Figure_S2/Panel_B_values.csv` | 73 | imc/workflows/compute_imc_signal_to_noise_qc.R |
| S20A | [Figure_S20/Panel_A.R](Figure_S20/Panel_A.R) | `Figure_S20/Panel_A_values.csv` | 20,657 | imc/workflows/Dist_ridge/Dist_to_cell_types_density_heatmap.Rmd |
| S20B | [Figure_S20/Panel_B.R](Figure_S20/Panel_B.R) | `Figure_S20/Panel_B_values.csv` | 20,657 | imc/workflows/Dist_ridge/Dist_to_cell_types_density_heatmap.Rmd |
| S20C | [Figure_S20/Panel_C.R](Figure_S20/Panel_C.R) | `Figure_S20/Panel_C_values.csv` | 32 | imc/workflows/Dist_ridge/Dist_to_cell_types_density_heatmap.Rmd |
| S20D | [Figure_S20/Panel_D.R](Figure_S20/Panel_D.R) | `Figure_S20/Panel_D_values.csv` | 32 | imc/workflows/Dist_ridge/Dist_to_cell_types_density_heatmap.Rmd |
| S21B | [Figure_S21/Panel_B.R](Figure_S21/Panel_B.R) | `Figure_S21/Panel_B_values.csv` | 4 | mintflow/workflows/redo_main_mintflow_cached_other_ec_labels_fast.R |
| S21C | [Figure_S21/Panel_C.R](Figure_S21/Panel_C.R) | `Figure_S21/Panel_C_values.csv` | 326 | mintflow/workflows/redo_main_mintflow_cached_other_ec_labels_fast.R |
| S22B | [Figure_S22/Panel_B.R](Figure_S22/Panel_B.R) | `Figure_S22/Panel_B_values.csv` | 6 | imc/workflows/Freq/Freq_clust.Rmd; imc/workflows/Met_expr/Met_expr.Rmd |
| S22C | [Figure_S22/Panel_C.R](Figure_S22/Panel_C.R) | `Figure_S22/Panel_C_values.csv` | 189 | functional/additional_source_workflows/cfc5d9ba_analyze_seahorse_crossed_design.R |
| S23A | [Figure_S23/Panel_A.R](Figure_S23/Panel_A.R) | `Figure_S23/Panel_A_values.csv` | 147 | imc/workflows/Met_expr/Met_expr.Rmd; imc/workflows/20231212_FibClustering.Rmd |
| S23B | [Figure_S23/Panel_B.R](Figure_S23/Panel_B.R) | `Figure_S23/Panel_B_values.csv` | 16 | imc/workflows/Freq/Freq_clust.Rmd; imc/workflows/Met_expr/Met_expr.Rmd |
| S23C | [Figure_S23/Panel_C.R](Figure_S23/Panel_C.R) | `Figure_S23/Panel_C_values.csv` | 16 | imc/workflows/Freq/Freq_clust.Rmd; imc/workflows/Met_expr/Met_expr.Rmd |
| S24A | [Figure_S24/Panel_A.R](Figure_S24/Panel_A.R) | `Figure_S24/Panel_A_values.csv` | 112 | imc/workflows/Met_expr/Met_expr_CN.Rmd |
| S24B | [Figure_S24/Panel_B.R](Figure_S24/Panel_B.R) | `Figure_S24/Panel_B_values.csv` | 70 | imc/workflows/Met_expr/Met_expr_CN.Rmd |
| S24C | [Figure_S24/Panel_C.R](Figure_S24/Panel_C.R) | `Figure_S24/Panel_C_values.csv` | 56 | imc/workflows/Met_expr/Met_expr_CN.Rmd |
| S25unlettered | [Figure_S25/Panel_unlettered.R](Figure_S25/Panel_unlettered.R) | `Figure_S25/Panel_unlettered_values.csv` | 578 | imc/workflows/Interaction/Interaction_tests_plots_beautiful_plots.Rmd |
| S26A | [Figure_S26/Panel_A.R](Figure_S26/Panel_A.R) | `Figure_S26/Panel_A_values.csv` | 8 | imc/workflows/Regression_clinicaldata/Differences_binary_clinical_outcome_dotplots.Rmd |
| S26B | [Figure_S26/Panel_B.R](Figure_S26/Panel_B.R) | `Figure_S26/Panel_B_values.csv` | 5 | imc/workflows/Regression_clinicaldata/Differences_binary_clinical_outcome_dotplots.Rmd |
| S26C | [Figure_S26/Panel_C.R](Figure_S26/Panel_C.R) | `Figure_S26/Panel_C_values.csv` | 4 | imc/workflows/Regression_clinicaldata/Differences_binary_clinical_outcome_dotplots.Rmd |
| S27A | [Figure_S27/Panel_A.R](Figure_S27/Panel_A.R) | `Figure_S27/Panel_A_values.csv` | 7 | imc/workflows/first_cohort_mrss_common.R; imc/workflows/Regression_clinicaldata/Clinical_correlations_dotplots.Rmd |
| S27B | [Figure_S27/Panel_B.R](Figure_S27/Panel_B.R) | `Figure_S27/Panel_B_values.csv` | 7 | imc/workflows/first_cohort_mrss_common.R; imc/workflows/Regression_clinicaldata/Clinical_correlations_dotplots.Rmd |
| S27C | [Figure_S27/Panel_C.R](Figure_S27/Panel_C.R) | `Figure_S27/Panel_C_values.csv` | 7 | imc/workflows/first_cohort_mrss_common.R; imc/workflows/Regression_clinicaldata/Clinical_correlations_dotplots.Rmd |
| S27D | [Figure_S27/Panel_D.R](Figure_S27/Panel_D.R) | `Figure_S27/Panel_D_values.csv` | 7 | imc/workflows/first_cohort_mrss_common.R; imc/workflows/Regression_clinicaldata/Clinical_correlations_dotplots.Rmd |
| S27E | [Figure_S27/Panel_E.R](Figure_S27/Panel_E.R) | `Figure_S27/Panel_E_values.csv` | 7 | imc/workflows/first_cohort_mrss_common.R; imc/workflows/Regression_clinicaldata/Clinical_correlations_dotplots.Rmd |
| S27F | [Figure_S27/Panel_F.R](Figure_S27/Panel_F.R) | `Figure_S27/Panel_F_values.csv` | 7 | imc/workflows/first_cohort_mrss_common.R; imc/workflows/Regression_clinicaldata/Clinical_correlations_dotplots.Rmd |
| S27G | [Figure_S27/Panel_G.R](Figure_S27/Panel_G.R) | `Figure_S27/Panel_G_values.csv` | 7 | imc/workflows/first_cohort_mrss_common.R; imc/workflows/Regression_clinicaldata/Clinical_correlations_dotplots.Rmd |
| S28A | [Figure_S28/Panel_A.R](Figure_S28/Panel_A.R) | `Figure_S28/Panel_A_values.csv` | 9 | imc/workflows/first_cohort_mrss_common.R; imc/workflows/Regression_clinicaldata/Clinical_correlations_dotplots.Rmd |
| S28B | [Figure_S28/Panel_B.R](Figure_S28/Panel_B.R) | `Figure_S28/Panel_B_values.csv` | 9 | imc/workflows/first_cohort_mrss_common.R; imc/workflows/Regression_clinicaldata/Clinical_correlations_dotplots.Rmd |
| S28C | [Figure_S28/Panel_C.R](Figure_S28/Panel_C.R) | `Figure_S28/Panel_C_values.csv` | 9 | imc/workflows/first_cohort_mrss_common.R; imc/workflows/Regression_clinicaldata/Clinical_correlations_dotplots.Rmd |
| S29A | [Figure_S29/Panel_A.R](Figure_S29/Panel_A.R) | `Figure_S29/Panel_A_values.csv` | 8 | imc/workflows/first_cohort_mrss_common.R; imc/workflows/Regression_clinicaldata/Clinical_correlations_dotplots.Rmd |
| S29B | [Figure_S29/Panel_B.R](Figure_S29/Panel_B.R) | `Figure_S29/Panel_B_values.csv` | 5 | imc/workflows/first_cohort_mrss_common.R; imc/workflows/Regression_clinicaldata/Clinical_correlations_dotplots.Rmd |
| S29C | [Figure_S29/Panel_C.R](Figure_S29/Panel_C.R) | `Figure_S29/Panel_C_values.csv` | 8 | imc/workflows/first_cohort_mrss_common.R; imc/workflows/Regression_clinicaldata/Clinical_correlations_dotplots.Rmd |
| S29D | [Figure_S29/Panel_D.R](Figure_S29/Panel_D.R) | `Figure_S29/Panel_D_values.csv` | 8 | imc/workflows/first_cohort_mrss_common.R; imc/workflows/Regression_clinicaldata/Clinical_correlations_dotplots.Rmd |
| S29E | [Figure_S29/Panel_E.R](Figure_S29/Panel_E.R) | `Figure_S29/Panel_E_values.csv` | 4 | imc/workflows/first_cohort_mrss_common.R; imc/workflows/Regression_clinicaldata/Clinical_correlations_dotplots.Rmd |
| S3A | [Figure_S3/Panel_A.R](Figure_S3/Panel_A.R) | `Figure_S3/Panel_A_values.csv` | 80 | imc/workflows/Freq/Freq_clust.Rmd; imc/workflows/Met_expr/Met_expr.Rmd |
| S3B | [Figure_S3/Panel_B.R](Figure_S3/Panel_B.R) | `Figure_S3/Panel_B_values.csv` | 80 | imc/workflows/Freq/Freq_clust.Rmd; imc/workflows/Met_expr/Met_expr.Rmd |
| S30A | [Figure_S30/Panel_A.R](Figure_S30/Panel_A.R) | `Figure_S30/Panel_A_values.csv` | 8 | imc/workflows/first_cohort_mrss_common.R; imc/workflows/Regression_clinicaldata/Clinical_correlations_dotplots.Rmd |
| S30B | [Figure_S30/Panel_B.R](Figure_S30/Panel_B.R) | `Figure_S30/Panel_B_values.csv` | 5 | imc/workflows/first_cohort_mrss_common.R; imc/workflows/Regression_clinicaldata/Clinical_correlations_dotplots.Rmd |
| S30C | [Figure_S30/Panel_C.R](Figure_S30/Panel_C.R) | `Figure_S30/Panel_C_values.csv` | 8 | imc/workflows/first_cohort_mrss_common.R; imc/workflows/Regression_clinicaldata/Clinical_correlations_dotplots.Rmd |
| S30D | [Figure_S30/Panel_D.R](Figure_S30/Panel_D.R) | `Figure_S30/Panel_D_values.csv` | 8 | imc/workflows/first_cohort_mrss_common.R; imc/workflows/Regression_clinicaldata/Clinical_correlations_dotplots.Rmd |
| S30E | [Figure_S30/Panel_E.R](Figure_S30/Panel_E.R) | `Figure_S30/Panel_E_values.csv` | 4 | imc/workflows/first_cohort_mrss_common.R; imc/workflows/Regression_clinicaldata/Clinical_correlations_dotplots.Rmd |
| S4A | [Figure_S4/Panel_A.R](Figure_S4/Panel_A.R) | `Figure_S4/Panel_A_values.csv` | 45,340 | imc/workflows/20231212_FibClustering.Rmd |
| S4B | [Figure_S4/Panel_B.R](Figure_S4/Panel_B.R) | `Figure_S4/Panel_B_values.csv` | 45,340 | imc/workflows/20231212_FibClustering.Rmd |
| S4C | [Figure_S4/Panel_C.R](Figure_S4/Panel_C.R) | `Figure_S4/Panel_C_values.csv` | 45,340 | imc/workflows/20231212_FibClustering.Rmd |
| S4D | [Figure_S4/Panel_D.R](Figure_S4/Panel_D.R) | `Figure_S4/Panel_D_values.csv` | 45,340 | imc/workflows/20231212_FibClustering.Rmd |
| S5A | [Figure_S5/Panel_A.R](Figure_S5/Panel_A.R) | `Figure_S5/Panel_A_values.csv` | 16 | imc/workflows/compute_cell_density_reviewer_qc.R |
| S5B | [Figure_S5/Panel_B.R](Figure_S5/Panel_B.R) | `Figure_S5/Panel_B_values.csv` | 16 | imc/workflows/compute_cell_density_reviewer_qc.R |
| S6A | [Figure_S6/Panel_A.R](Figure_S6/Panel_A.R) | `Figure_S6/Panel_A_values.csv` | 20,657 | imc/workflows/20231212_FibClustering.Rmd |
| S6B | [Figure_S6/Panel_B.R](Figure_S6/Panel_B.R) | `Figure_S6/Panel_B_values.csv` | 105 | imc/workflows/Met_expr/Met_expr.Rmd; imc/workflows/20231212_FibClustering.Rmd |
| S6C | [Figure_S6/Panel_C.R](Figure_S6/Panel_C.R) | `Figure_S6/Panel_C_values.csv` | 3,288 | imc/workflows/20231212_FibClustering.Rmd; shared/workflows/IMC_heatmap_sankey_new.Rmd |
| S7A | [Figure_S7/Panel_A.R](Figure_S7/Panel_A.R) | `Figure_S7/Panel_A_values.csv` | 20,657 | imc/workflows/20231212_FibClustering.Rmd; shared/workflows/IMC_heatmap_sankey_new.Rmd |
| S7B | [Figure_S7/Panel_B.R](Figure_S7/Panel_B.R) | `Figure_S7/Panel_B_values.csv` | 16 | imc/workflows/Freq/Freq_clust.Rmd |
| S7C | [Figure_S7/Panel_C.R](Figure_S7/Panel_C.R) | `Figure_S7/Panel_C_values.csv` | 16 | imc/workflows/Freq/Freq_clust.Rmd; imc/workflows/Met_expr/Met_expr.Rmd |
| S7D | [Figure_S7/Panel_D.R](Figure_S7/Panel_D.R) | `Figure_S7/Panel_D_values.csv` | 16 | imc/workflows/Freq/Freq_clust.Rmd; imc/workflows/Met_expr/Met_expr.Rmd |
| S7E | [Figure_S7/Panel_E.R](Figure_S7/Panel_E.R) | `Figure_S7/Panel_E_values.csv` | 16 | imc/workflows/Freq/Freq_clust.Rmd; imc/workflows/Met_expr/Met_expr.Rmd |
| S7F | [Figure_S7/Panel_F.R](Figure_S7/Panel_F.R) | `Figure_S7/Panel_F_values.csv` | 16 | imc/workflows/Freq/Freq_clust.Rmd; imc/workflows/Met_expr/Met_expr.Rmd |
| S7G | [Figure_S7/Panel_G.R](Figure_S7/Panel_G.R) | `Figure_S7/Panel_G_values.csv` | 16 | imc/workflows/Freq/Freq_clust.Rmd; imc/workflows/Met_expr/Met_expr.Rmd |
| S7H | [Figure_S7/Panel_H.R](Figure_S7/Panel_H.R) | `Figure_S7/Panel_H_values.csv` | 16 | imc/workflows/Freq/Freq_clust.Rmd; imc/workflows/Met_expr/Met_expr.Rmd |
| S7I | [Figure_S7/Panel_I.R](Figure_S7/Panel_I.R) | `Figure_S7/Panel_I_values.csv` | 16 | imc/workflows/Freq/Freq_clust.Rmd; imc/workflows/Met_expr/Met_expr.Rmd |
| S7J | [Figure_S7/Panel_J.R](Figure_S7/Panel_J.R) | `Figure_S7/Panel_J_values.csv` | 16 | imc/workflows/Freq/Freq_clust.Rmd; imc/workflows/Met_expr/Met_expr.Rmd |
| S8A | [Figure_S8/Panel_A.R](Figure_S8/Panel_A.R) | `Figure_S8/Panel_A_values.csv` | 19,183 | imc/workflows/20231212_FibClustering.Rmd; shared/workflows/IMC_heatmap_sankey_new.Rmd |
| S8B | [Figure_S8/Panel_B.R](Figure_S8/Panel_B.R) | `Figure_S8/Panel_B_values.csv` | 128 | shared/workflows/IMC_heatmap_sankey_new.Rmd |
| S9A | [Figure_S9/Panel_A.R](Figure_S9/Panel_A.R) | `Figure_S9/Panel_A_values.csv` | 42 | shared/workflows/six_donor_imc_only_figure_panels.Rmd |
| S9B | [Figure_S9/Panel_B.R](Figure_S9/Panel_B.R) | `Figure_S9/Panel_B_values.csv` | 24 | shared/workflows/six_donor_imc_only_figure_panels.Rmd |
| S9C | [Figure_S9/Panel_C.R](Figure_S9/Panel_C.R) | `Figure_S9/Panel_C_values.csv` | 950 | shared/workflows/IMC_protein_and_scMetabolism_scores_combined.Rmd |
| S9D | [Figure_S9/Panel_D.R](Figure_S9/Panel_D.R) | `Figure_S9/Panel_D_values.csv` | 21,812 | shared/reference_projection.R |
| S9E | [Figure_S9/Panel_E.R](Figure_S9/Panel_E.R) | `Figure_S9/Panel_E_values.csv` | 3,105 | imc/workflows/20231212_FibClustering.Rmd; shared/workflows/IMC_heatmap_sankey_new.Rmd |
| S9F | [Figure_S9/Panel_F.R](Figure_S9/Panel_F.R) | `Figure_S9/Panel_F_values.csv` | 128 | shared/workflows/IMC_heatmap_sankey_new.Rmd |

Image/schematic panel scripts use `ASSET_DIR/Figure_<id>/Panel_<id>_image.png`. They export original pixels instead of synthesizing images. These panels are listed in the deposited `Panel_index.csv` with no numerical table. Mixed microscopy/quantification panels also retain original microscopy assets for manuscript composition.
