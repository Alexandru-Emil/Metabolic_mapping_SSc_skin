script_dir <- {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)) else getwd()
}

stages <- c(
  "01_qc_batch_correction.R",
  "02_deg_edger.R",
  "03_functional_analysis.R",
  "05_method_threshold_sensitivity.R",
  "06_combatseq_deg_go_reports.R",
  "07_symbol_interactive_reports.R",
  "08_threshold_go_gsea_reports.R",
  "04_report.R"
)

for (stage in stages) {
  message("\n--- Running ", stage, " ---")
  stage_path <- normalizePath(file.path(script_dir, stage), winslash = "/", mustWork = TRUE)
  status <- system2("Rscript", c(shQuote(stage_path), commandArgs(trailingOnly = TRUE)))
  message("Stage status: ", status)
  if (!is.null(status) && status != 0) {
    stop("Stage failed: ", stage, call. = FALSE)
  }
}

message("\nAll edgeR PCS-batch-corrected RNA-seq stages complete.")
