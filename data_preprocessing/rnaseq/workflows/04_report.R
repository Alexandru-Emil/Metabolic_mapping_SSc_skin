# Scientific workflow derived from Revision/RNAseq PCS met inhibitors/results/scripts/04_report.R
# Run from the repository root with METABOLIC_INPUT_DIR set to the external project data.
source("data_preprocessing/common/config.R")
script_dir <- {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)) else getwd()
}
source(file.path(script_dir, "00_config.R"))

stage_name <- "04_report"
stage_rds <- file.path(cache_dir, paste0(stage_name, ".rds"))
stage1_rds <- file.path(cache_dir, "01_qc_batch_correction.rds")
stage2_rds <- file.path(cache_dir, "02_deg_edger.rds")
stage3_rds <- file.path(cache_dir, "03_functional_analysis.rds")
stage5_rds <- file.path(cache_dir, "05_method_threshold_sensitivity.rds")
deps <- c(stage1_rds, stage2_rds, stage3_rds, stage5_rds, file.path(script_dir, "00_config.R"), file.path(script_dir, "04_report.R"))

if (stage_is_fresh(stage_rds, deps)) {
  message("Using cached ", stage_name, ": ", stage_rds)
  quit(save = "no", status = 0)
}

message("Stage 4: workbook and HTML technical report")
s1 <- read_stage(stage1_rds)
s2 <- read_stage(stage2_rds)
s3 <- read_stage(stage3_rds)
s5 <- read_stage(stage5_rds)

wb <- createWorkbook()
addWorksheet(wb, "sample_metadata")
writeDataTable(wb, "sample_metadata", s1$sample_meta)
addWorksheet(wb, "lane_metadata")
writeDataTable(wb, "lane_metadata", s1$lane_meta)
addWorksheet(wb, "filtering_summary")
writeDataTable(wb, "filtering_summary", s1$filtering_summary)
addWorksheet(wb, "rounding_summary")
writeDataTable(wb, "rounding_summary", s1$rounding_summary)
addWorksheet(wb, "edgeR_summary")
writeDataTable(wb, "edgeR_summary", s2$deg_summary)
for (nm in names(s2$edgeR_results)) {
  sheet <- paste0("edgeR_", substr(nm, 1, 24))
  addWorksheet(wb, sheet)
  writeDataTable(wb, sheet, s2$edgeR_results[[nm]])
  sig_sheet <- paste0("DEG_", substr(nm, 1, 25))
  addWorksheet(wb, sig_sheet)
  writeDataTable(wb, sig_sheet, s2$edgeR_results[[nm]][s2$edgeR_results[[nm]]$FDR < 0.05 & abs(s2$edgeR_results[[nm]]$logFC) >= 1, ])
  gsea_sheet <- paste0("GSEA_", substr(nm, 1, 24))
  addWorksheet(wb, gsea_sheet)
  writeDataTable(wb, gsea_sheet, s3$gsea_results[[nm]])
}
addWorksheet(wb, "method_threshold_key")
writeDataTable(wb, "method_threshold_key", s5$key_threshold_summary)
addWorksheet(wb, "method_status")
writeDataTable(wb, "method_status", s5$method_status)
addWorksheet(wb, "top_multimethod")
writeDataTable(wb, "top_multimethod", s5$top_multimethod_candidates)
addWorksheet(wb, "no_deg_diagnostics")
writeDataTable(wb, "no_deg_diagnostics", s5$primary_diagnostics)
addWorksheet(wb, "dispersion_power")
writeDataTable(wb, "dispersion_power", s5$dispersion_summary)
addWorksheet(wb, "variance_partition")
writeDataTable(wb, "variance_partition", s5$variance_partition_summary)
addWorksheet(wb, "possible_reasons")
writeDataTable(wb, "possible_reasons", s5$possible_reasons)
saveWorkbook(wb, file.path(table_dir, "bulk_rnaseq_edgeR_PCS_batch_corrected_publication_results.xlsx"), overwrite = TRUE)

summary_obj <- list(
  source_notes = list(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    input_file = input_file,
    analysis_dir = analysis_dir,
    qc_stage = "01_qc_batch_correction.R: lane aggregation, TMM normalization, ComBat-seq PCS batch correction, QC figures",
    deg_stage = "02_deg_edger.R: edgeR QL GLM with design ~ PCS + condition",
    functional_stage = "03_functional_analysis.R: GO BP GSEA/ORA from edgeR model statistics",
    sensitivity_stage = "05_method_threshold_sensitivity.R: Shiny-app-style method and p/FDR/fold-change threshold sensitivity",
    report_stage = "04_report.R: workbook and HTML report",
    batch_correction_method = s1$batch_correction_method,
    caveats = c(
      "PCS was modeled as a fixed batch/block covariate in edgeR DEG testing.",
      "ComBat-seq PCS-corrected counts/logCPM are used for downstream visualization and expression summaries, not as edgeR model input.",
      "The source matrix contains fractional estimated counts; rounded lane-summed counts are used for edgeR and ComBat-seq."
    )
  ),
  filtering_summary = s1$filtering_summary,
  rounding_summary = s1$rounding_summary,
  sample_meta = s1$sample_meta,
  deg_summary = s2$deg_summary,
  top_genes_summary = s2$top_genes_summary,
  top_gsea_summary = s3$top_gsea_summary,
  method_threshold_key_summary = s5$key_threshold_summary,
  top_multimethod_candidates = s5$top_multimethod_candidates,
  primary_diagnostics = s5$primary_diagnostics,
  dispersion_summary = s5$dispersion_summary,
  variance_partition_summary = s5$variance_partition_summary,
  possible_reasons = s5$possible_reasons,
  method_status = s5$method_status
)
saveRDS(summary_obj, file.path(object_dir, "edgeR_batch_corrected_analysis_summary.rds"))

rmd_path <- file.path(report_dir, "bulk_rnaseq_edgeR_PCS_batch_corrected_report.Rmd")
report_lines <- c(
  "---",
  "title: \"Bulk RNA-seq edgeR analysis with PCS batch correction and DEG sensitivity\"",
  "output:",
  "  html_document:",
  "    toc: true",
  "    toc_depth: 3",
  "    number_sections: false",
  "    theme: cosmo",
  "    df_print: paged",
  "params:",
  "  out_dir: \"\"",
  "---",
  "",
  "```{r setup, include=FALSE}",
  "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
  "out_dir <- params$out_dir",
  "summary <- readRDS(file.path(out_dir, 'objects', 'edgeR_batch_corrected_analysis_summary.rds'))",
  "fmt_int <- function(x) format(x, big.mark = ',', scientific = FALSE, trim = TRUE)",
  "fig <- function(name) knitr::include_graphics(file.path(out_dir, 'figures', name))",
  "deg_summary <- summary$deg_summary",
  "top_genes <- summary$top_genes_summary",
  "top_gsea <- summary$top_gsea_summary",
  "method_thresholds <- summary$method_threshold_key_summary",
  "top_multi <- summary$top_multimethod_candidates",
  "primary_diag <- summary$primary_diagnostics",
  "reasons <- summary$possible_reasons",
  "method_status <- summary$method_status",
  "filtering <- summary$filtering_summary",
  "sample_meta <- summary$sample_meta",
  "```",
  "",
  "## Technical summary",
  "",
  "The previous analysis accounted for PCS identity by including PCS as a paired fixed effect in the model, but it did not create a batch-corrected expression matrix. This refactored pipeline treats each PCS block as a batch in two complementary ways: PCS is included directly in the edgeR quasi-likelihood GLM for valid differential-expression inference, and PCS-corrected expression matrices are exported for downstream visualization and expression summaries.",
  "",
  "The source matrix contains fractional estimated counts, so technical lanes were first summed to biological PCS-condition samples and then rounded for the count-based edgeR/ComBat-seq stages. `r fmt_int(filtering$value[filtering$metric == 'genes_after_filterByExpr'])` genes were retained for testing across `r fmt_int(filtering$value[filtering$metric == 'biological_samples'])` biological profiles. High-confidence DEG calls use FDR < 0.05 and |log2FC| >= 1.",
  "",
  "```{r summary-table}",
  "show_summary <- deg_summary[, c('contrast_label', 'genes_tested', 'fdr_0_05', 'up_fdr_0_05_abslog2fc_1', 'down_fdr_0_05_abslog2fc_1', 'deg_fdr_0_05_abslog2fc_1')]",
  "colnames(show_summary) <- c('Contrast', 'Genes tested', 'FDR < 0.05', 'Up DEG', 'Down DEG', 'Total DEG')",
  "knitr::kable(show_summary, digits = 3)",
  "```",
  "",
  "## Shiny-app-style method and threshold sensitivity",
  "",
  "The accompanying Shiny app supports edgeR exact/classical workflows, DESeq2 workflows, condition-preserving ComBat-seq batch adjustment, and user-selected p-value/FDR plus fold-change filters. This report mirrors those ideas in a reproducible cached stage while keeping the primary manuscript inference as edgeR QL with PCS in the design. Methods that operate directly on ComBat-seq corrected counts or corrected logCPM are interpreted as exploratory candidate screens.",
  "",
  "```{r method-status}",
  "ms <- method_status",
  "colnames(ms) <- c('Method', 'Status', 'Detail')",
  "knitr::kable(ms)",
  "```",
  "",
  "```{r threshold-heatmap, out.width='100%'}",
  "fig('deg_method_threshold_sensitivity_heatmap.png')",
  "```",
  "",
  "The table below summarizes the key thresholds most similar to the Shiny app controls. `FC >= 1.5` means |log2FC| >= log2(1.5), while `FC >= 2` is equivalent to |log2FC| >= 1.",
  "",
  "```{r key-thresholds}",
  "kt <- method_thresholds[, c('method', 'contrast_label', 'threshold_metric', 'alpha', 'fold_change', 'up', 'down', 'total_deg', 'min_pvalue', 'min_fdr')]",
  "colnames(kt) <- c('Method', 'Contrast', 'Metric', 'Alpha', 'Fold-change', 'Up', 'Down', 'Total DEG', 'Minimum P', 'Minimum FDR')",
  "knitr::kable(kt, digits = 4)",
  "```",
  "",
  "Top candidate tables are exported in full. The preview below emphasizes genes passing at least one relaxed or publication-like threshold where available, followed by the strongest ranked genes.",
  "",
  "```{r top-multimethod}",
  "tm <- top_multi[, c('method', 'contrast_label', 'selection_tier', 'gene_label', 'GENENAME', 'logFC', 'PValue', 'FDR')]",
  "colnames(tm) <- c('Method', 'Contrast', 'Selection tier', 'Gene', 'Gene name', 'log2FC', 'P value', 'FDR')",
  "knitr::kable(head(tm, 50), digits = 4)",
  "```",
  "",
  "## PCS correction changes the exploratory expression space, not the count-model contract",
  "",
  "The PCA compares TMM logCPM before and after PCS correction. The corrected panel uses ComBat-seq count correction with condition preserved, then TMM logCPM for visualization. The edgeR DEG model still uses rounded lane-summed counts with PCS in the design, which avoids fitting a negative-binomial count model to transformed values.",
  "",
  "```{r pca, out.width='100%'}",
  "fig('qc_pca_before_after_PCS_batch_correction.png')",
  "```",
  "",
  "Library sizes and sample correlations provide the basic QC envelope for the edgeR input and the PCS-corrected downstream matrix.",
  "",
  "```{r lib, out.width='100%'}",
  "fig('qc_library_sizes_edgeR_inputs.png')",
  "```",
  "",
  "```{r cor, out.width='100%'}",
  "fig('qc_sample_correlation_after_PCS_batch_correction.png')",
  "```",
  "",
  "## edgeR estimates treatment effects while adjusting for PCS batch",
  "",
  "The primary model is an edgeR quasi-likelihood GLM with TMM normalization and design `~ PCS + condition`. The treatment coefficients test 2DG + OX vs control and LDHA + OX vs control while controlling for PCS block.",
  "",
  "```{r deg-bar, out.width='100%'}",
  "fig('deg_edgeR_summary_bar.png')",
  "```",
  "",
  "### 2DG + OX versus control",
  "",
  "```{r volcano-2dg, out.width='100%'}",
  "fig('deg_edgeR_volcano_2DG_vs_C.png')",
  "```",
  "",
  "```{r ma-2dg, out.width='100%'}",
  "fig('deg_edgeR_ma_2DG_vs_C.png')",
  "```",
  "",
  "### LDHA + OX versus control",
  "",
  "```{r volcano-ldha, out.width='100%'}",
  "fig('deg_edgeR_volcano_LDHA_vs_C.png')",
  "```",
  "",
  "```{r ma-ldha, out.width='100%'}",
  "fig('deg_edgeR_ma_LDHA_vs_C.png')",
  "```",
  "",
  "## PCS-corrected expression supports the heatmap shortlist",
  "",
  "The heatmap uses PCS-corrected logCPM values for top edgeR-ranked genes across both contrasts. Exact statistics should be taken from the edgeR all-gene tables.",
  "",
  "```{r heatmap, out.width='100%'}",
  "fig('deg_edgeR_top_genes_PCS_corrected_heatmap.png')",
  "```",
  "",
  "```{r top-genes}",
  "tg <- top_genes[, c('contrast_label', 'gene_label', 'GENENAME', 'logFC', 'logCPM', 'FDR', 'direction')]",
  "colnames(tg) <- c('Contrast', 'Gene', 'Gene name', 'log2FC', 'logCPM', 'FDR', 'Class')",
  "knitr::kable(tg, digits = 4)",
  "```",
  "",
  "## GO Biological Process analysis uses edgeR model statistics",
  "",
  "GSEA was run on signed square-root edgeR QL F-statistics, so pathway ranking reflects treatment effects after PCS adjustment in the count model. Positive NES values indicate enrichment toward treatment-up genes; negative NES values indicate enrichment toward control-up genes.",
  "",
  "```{r gsea-2dg, out.width='100%'}",
  "fig('pathway_edgeR_go_bp_gsea_2DG_vs_C.png')",
  "```",
  "",
  "```{r gsea-ldha, out.width='100%'}",
  "fig('pathway_edgeR_go_bp_gsea_LDHA_vs_C.png')",
  "```",
  "",
  "```{r top-gsea}",
  "if (nrow(top_gsea) > 0) {",
  "  gsea_show <- top_gsea[, c('contrast_label', 'Description', 'setSize', 'NES', 'p.adjust', 'qvalue')]",
  "  colnames(gsea_show) <- c('Contrast', 'GO Biological Process', 'Set size', 'NES', 'FDR', 'q value')",
  "  knitr::kable(gsea_show, digits = 4)",
  "} else {",
  "  cat('No GO GSEA results were returned.')",
  "}",
  "```",
  "",
  "## Why high-confidence DEGs may be sparse",
  "",
  "The absence or scarcity of high-confidence single-gene calls is not, by itself, evidence that the experiment failed. In this design, the most defensible DEG model is conservative: it fits PCS block effects, estimates gene-wise negative-binomial dispersion, and corrects thousands of simultaneous tests. The diagnostics below separate statistical power limits from exploratory candidate signals.",
  "",
  "```{r primary-diagnostics}",
  "pd <- primary_diag[, c('contrast_label', 'genes_tested', 'min_pvalue', 'min_fdr', 'raw_p_lt_0_05', 'raw_p_lt_0_05_and_fc_ge_1_5', 'raw_p_lt_0_05_and_fc_ge_2', 'top_gene', 'top_gene_logFC', 'top_gene_PValue', 'top_gene_FDR')]",
  "colnames(pd) <- c('Contrast', 'Genes tested', 'Minimum P', 'Minimum FDR', 'Raw P<0.05', 'Raw P<0.05 and FC>=1.5', 'Raw P<0.05 and FC>=2', 'Top gene', 'Top gene log2FC', 'Top gene P', 'Top gene FDR')",
  "knitr::kable(pd, digits = 4)",
  "```",
  "",
  "```{r reasons-table}",
  "rs <- reasons[, c('reason', 'evidence', 'interpretation')]",
  "colnames(rs) <- c('Reason', 'Evidence in this dataset', 'Interpretation')",
  "knitr::kable(rs)",
  "```",
  "",
  "## Scope, data, and metric definitions",
  "",
  "- Biological comparisons: 2DG + OX vs control and LDHA + OX vs control.",
  "- Batch/block variable: PCS identity, modeled in edgeR and corrected in exported downstream expression matrices.",
  "- Count input for edgeR: technical-lane-summed counts rounded to integers because the source matrix contains fractional estimated counts.",
  "- Batch-corrected downstream matrix: ComBat-seq PCS-corrected counts and derived TMM logCPM, with condition preserved.",
  "- DEG threshold: FDR < 0.05 and |log2FC| >= 1 for publication shortlist; all tested genes are exported.",
  "- Sensitivity thresholds: Shiny-style filters include FDR or raw P value cutoffs at 0.05, 0.10, and 0.20 with fold-change cutoffs of 1, 1.2, 1.5, and 2.",
  "",
  "```{r sample-table}",
  "sm <- sample_meta[, c('sample_id', 'pcs', 'condition_label', 'n_technical_lanes', 'library_size')]",
  "colnames(sm) <- c('Sample', 'PCS block', 'Condition', 'Technical lanes', 'Summed library size')",
  "knitr::kable(sm)",
  "```",
  "",
  "## Methodology",
  "",
  "1. `01_qc_batch_correction.R` imports the Galaxy count matrix, parses lane metadata, sums technical lanes, annotates genes, filters low-expression genes, performs TMM normalization, and creates PCS-corrected downstream matrices.",
  "2. `02_deg_edger.R` fits edgeR quasi-likelihood GLMs using `~ PCS + condition` and exports contrast-specific all-gene and DEG tables.",
  "3. `03_functional_analysis.R` runs GO BP GSEA on edgeR signed square-root QL F-statistics and ORA on up/down DEG sets when enough genes are available.",
  "4. `05_method_threshold_sensitivity.R` runs Shiny-app-style method sensitivity and p/FDR/fold-change threshold grids, then writes candidate and no-DEG diagnostic tables.",
  "5. `04_report.R` builds the Excel workbook and this HTML technical report.",
  "6. `run_all_edgeR_batch_corrected.R` orchestrates the stages. Each stage caches an RDS object and skips recomputation when inputs and scripts are unchanged; use `--force` to rebuild.",
  "",
  "## Limitations, uncertainty, and robustness checks",
  "",
  "- Batch-corrected values are exported for downstream visualization/expression summaries, but transformed or adjusted matrices are not used as edgeR model input.",
  "- The multi-method threshold grid is a sensitivity analysis. FDR-controlled edgeR QL results with PCS in the model should remain the primary manuscript statistics.",
  "- ComBat-seq correction treats PCS as batch and condition as the group to preserve. With five PCS blocks and one sample per condition per PCS, the correction is useful for visualization but should not be over-interpreted as removing every paired-sample effect.",
  "- The source matrix contains fractional estimated counts; rounding is necessary for edgeR/ComBat-seq count workflows and is documented in `rounding_summary_for_edgeR.csv`.",
  "- No upstream FASTQ/BAM QC, alignment metrics, strandedness checks, or gene-body coverage metrics were supplied.",
  "",
  "## Recommended next steps",
  "",
  "1. Use the edgeR all-gene tables as the primary DEG source for manuscript statistics.",
  "2. Use the PCS-corrected logCPM matrix for heatmaps, PCA-style exploratory figures, and focused expression panels.",
  "3. Add upstream sequencing/alignment QC if available before final submission.",
  "4. Review the GO terms and leading-edge genes manually to collapse redundant pathway labels into a concise biological narrative.",
  "",
  "## Further questions",
  "",
  "- Are there non-PCS covariates such as library-prep batch or RNA quality that should be included in the edgeR design?",
  "- Should batch-corrected expression panels show all genes of interest or only edgeR-ranked genes?",
  "- Are curated metabolic or perturbation-specific gene sets available to supplement GO?"
)
writeLines(report_lines, rmd_path)

Sys.setenv(RSTUDIO_PANDOC = project_path("external/tools"))
rmarkdown::render(
  input = rmd_path,
  output_file = "bulk_rnaseq_edgeR_PCS_batch_corrected_report.html",
  output_dir = report_dir,
  params = list(out_dir = analysis_dir),
  quiet = TRUE,
  envir = new.env(parent = globalenv())
)

saveRDS(list(timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), stage_rds)
write_stage_stamp(stage_name)
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_04_report.txt"))
message("Stage 4 complete: ", file.path(report_dir, "bulk_rnaseq_edgeR_PCS_batch_corrected_report.html"))
