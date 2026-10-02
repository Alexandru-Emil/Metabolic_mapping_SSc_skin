# Scientific workflow derived from Revision/RNAseq PCS met inhibitors/results/scripts/05_method_threshold_sensitivity.R
# Run from the repository root with METABOLIC_INPUT_DIR set to the external project data.
source("data_preprocessing/common/config.R")
script_dir <- {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)) else getwd()
}
source(file.path(script_dir, "00_config.R"))

stage_name <- "05_method_threshold_sensitivity"
stage_rds <- file.path(cache_dir, paste0(stage_name, ".rds"))
stage1_rds <- file.path(cache_dir, "01_qc_batch_correction.rds")
stage2_rds <- file.path(cache_dir, "02_deg_edger.rds")
deps <- c(
  stage1_rds,
  stage2_rds,
  file.path(script_dir, "00_config.R"),
  file.path(script_dir, "05_method_threshold_sensitivity.R")
)

if (stage_is_fresh(stage_rds, deps)) {
  message("Using cached ", stage_name, ": ", stage_rds)
  quit(save = "no", status = 0)
}

message("Stage 5: Shiny-app-style DEG method and threshold sensitivity analysis")
s1 <- read_stage(stage1_rds)
s2 <- read_stage(stage2_rds)

sample_meta <- s1$sample_meta
sample_meta$condition <- relevel(sample_meta$condition, ref = "C")
sample_meta$pcs <- droplevels(sample_meta$pcs)
sample_meta <- sample_meta[sample_meta$sample_id, , drop = FALSE]
annotation <- s1$filtered_annotation
annotation_small <- annotation[, intersect(c("gene_id", "ensembl_id", "SYMBOL", "GENENAME", "ENTREZID", "gene_label"), names(annotation)), drop = FALSE]

method_status <- data.frame(method = character(), status = character(), detail = character(), stringsAsFactors = FALSE)
add_status <- function(method, status, detail = "") {
  method_status <<- rbind(method_status, data.frame(method = method, status = status, detail = detail, stringsAsFactors = FALSE))
}

standardize_result <- function(df, method, contrast, contrast_label, logfc_col, p_col, fdr_col,
                               stat_col = NULL, expr_col = NULL, model_note = "") {
  if (is.null(df) || nrow(df) == 0) return(data.frame())
  if (!"gene_id" %in% names(df)) df$gene_id <- rownames(df)
  stat <- if (!is.null(stat_col) && stat_col %in% names(df)) df[[stat_col]] else rep(NA_real_, nrow(df))
  expr <- if (!is.null(expr_col) && expr_col %in% names(df)) df[[expr_col]] else rep(NA_real_, nrow(df))
  out <- data.frame(
    method = method,
    contrast = contrast,
    contrast_label = contrast_label,
    gene_id = as.character(df$gene_id),
    logFC = as.numeric(df[[logfc_col]]),
    statistic = as.numeric(stat),
    PValue = as.numeric(df[[p_col]]),
    FDR = as.numeric(df[[fdr_col]]),
    average_expression = as.numeric(expr),
    model_note = model_note,
    stringsAsFactors = FALSE
  )
  out <- merge(out, annotation_small, by = "gene_id", all.x = TRUE, sort = FALSE)
  out <- out[order(out$PValue, out$FDR, -abs(out$logFC), na.last = TRUE), ]
  out$rank_metric_signed_logP <- sign(out$logFC) * -log10(pmax(out$PValue, .Machine$double.xmin))
  out$gene_label <- ifelse(is.na(out$gene_label) | out$gene_label == "", out$gene_id, out$gene_label)
  rownames(out) <- NULL
  out
}

run_safely <- function(method, expr) {
  tryCatch({
    out <- force(expr)
    add_status(method, "ok", paste0(nrow(out), " standardized gene rows"))
    out
  }, error = function(e) {
    msg <- conditionMessage(e)
    add_status(method, "failed", msg)
    writeLines(msg, file.path(log_dir, paste0("method_sensitivity_error_", gsub("[^A-Za-z0-9]+", "_", method), ".txt")))
    data.frame()
  })
}

method_results <- list()

method_results[["edgeR_QL_PCS"]] <- run_safely("edgeR_QL_PCS", {
  do.call(rbind, lapply(names(s2$edgeR_results), function(nm) {
    standardize_result(
      s2$edgeR_results[[nm]],
      method = "edgeR_QL_PCS",
      contrast = nm,
      contrast_label = contrast_specs[[nm]]$label,
      logfc_col = "logFC",
      p_col = "PValue",
      fdr_col = "FDR",
      stat_col = "F",
      expr_col = "logCPM",
      model_note = "Primary publication model: edgeR quasi-likelihood GLM, TMM normalization, design ~ PCS + condition."
    )
  }))
})

method_results[["edgeR_LRT_PCS"]] <- run_safely("edgeR_LRT_PCS", {
  y_lrt <- s2$y
  fit_lrt <- glmFit(y_lrt, s2$design)
  do.call(rbind, lapply(names(contrast_specs), function(nm) {
    lrt <- glmLRT(fit_lrt, coef = contrast_specs[[nm]]$coef)
    tt <- as.data.frame(topTags(lrt, n = Inf, sort.by = "PValue"))
    tt$gene_id <- rownames(tt)
    standardize_result(
      tt,
      method = "edgeR_LRT_PCS",
      contrast = nm,
      contrast_label = contrast_specs[[nm]]$label,
      logfc_col = "logFC",
      p_col = "PValue",
      fdr_col = "FDR",
      stat_col = "LR",
      expr_col = "logCPM",
      model_note = "Sensitivity model: edgeR likelihood-ratio GLM, TMM normalization, design ~ PCS + condition."
    )
  }))
})

method_results[["edgeR_exact_ComBatSeq"]] <- run_safely("edgeR_exact_ComBatSeq", {
  bc_counts <- round(s1$batch_corrected_pseudo_counts[, sample_meta$sample_id, drop = FALSE])
  bc_counts[bc_counts < 0 | is.na(bc_counts)] <- 0
  storage.mode(bc_counts) <- "integer"
  group_all <- droplevels(sample_meta$condition)
  do.call(rbind, lapply(names(contrast_specs), function(nm) {
    cond <- contrast_specs[[nm]]$condition
    idx <- group_all %in% c("C", cond)
    y_pair <- DGEList(counts = bc_counts[, idx, drop = FALSE], group = droplevels(group_all[idx]))
    y_pair <- calcNormFactors(y_pair, method = "TMM")
    y_pair <- estimateDisp(y_pair, robust = TRUE)
    et <- exactTest(y_pair, pair = c("C", cond))
    tt <- as.data.frame(topTags(et, n = Inf, sort.by = "PValue"))
    tt$gene_id <- rownames(tt)
    standardize_result(
      tt,
      method = "edgeR_exact_ComBatSeq",
      contrast = nm,
      contrast_label = contrast_specs[[nm]]$label,
      logfc_col = "logFC",
      p_col = "PValue",
      fdr_col = "FDR",
      expr_col = "logCPM",
      model_note = "Exploratory Shiny-style exact test on ComBat-seq PCS-corrected counts; condition-only two-group comparison."
    )
  }))
})

if (requireNamespace("DESeq2", quietly = TRUE)) {
  method_results[["DESeq2_PCS"]] <- run_safely("DESeq2_PCS", {
    counts <- round(s1$filtered_counts[, sample_meta$sample_id, drop = FALSE])
    counts[counts < 0 | is.na(counts)] <- 0
    storage.mode(counts) <- "integer"
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = counts,
      colData = sample_meta,
      design = ~ pcs + condition
    )
    dds <- dds[rowSums(DESeq2::counts(dds)) > 0, ]
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    do.call(rbind, lapply(names(contrast_specs), function(nm) {
      cond <- contrast_specs[[nm]]$condition
      res <- as.data.frame(DESeq2::results(dds, contrast = c("condition", cond, "C"), alpha = 0.05))
      res$gene_id <- rownames(res)
      standardize_result(
        res,
        method = "DESeq2_PCS",
        contrast = nm,
        contrast_label = contrast_specs[[nm]]$label,
        logfc_col = "log2FoldChange",
        p_col = "pvalue",
        fdr_col = "padj",
        stat_col = "stat",
        expr_col = "baseMean",
        model_note = "Sensitivity model: DESeq2 Wald test, design ~ PCS + condition."
      )
    }))
  })

  method_results[["DESeq2_ComBatSeq"]] <- run_safely("DESeq2_ComBatSeq", {
    bc_counts <- round(s1$batch_corrected_pseudo_counts[, sample_meta$sample_id, drop = FALSE])
    bc_counts[bc_counts < 0 | is.na(bc_counts)] <- 0
    storage.mode(bc_counts) <- "integer"
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = bc_counts,
      colData = sample_meta,
      design = ~ condition
    )
    dds <- dds[rowSums(DESeq2::counts(dds)) > 0, ]
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    do.call(rbind, lapply(names(contrast_specs), function(nm) {
      cond <- contrast_specs[[nm]]$condition
      res <- as.data.frame(DESeq2::results(dds, contrast = c("condition", cond, "C"), alpha = 0.05))
      res$gene_id <- rownames(res)
      standardize_result(
        res,
        method = "DESeq2_ComBatSeq",
        contrast = nm,
        contrast_label = contrast_specs[[nm]]$label,
        logfc_col = "log2FoldChange",
        p_col = "pvalue",
        fdr_col = "padj",
        stat_col = "stat",
        expr_col = "baseMean",
        model_note = "Exploratory Shiny-style DESeq2 Wald test on ComBat-seq PCS-corrected counts; design ~ condition."
      )
    }))
  })
} else {
  add_status("DESeq2_PCS", "skipped", "DESeq2 package is not installed.")
  add_status("DESeq2_ComBatSeq", "skipped", "DESeq2 package is not installed.")
}

method_results[["limma_trend_ComBatSeq_logCPM"]] <- run_safely("limma_trend_ComBatSeq_logCPM", {
  expr <- s1$batch_corrected_logcpm[, sample_meta$sample_id, drop = FALSE]
  design_condition <- model.matrix(~ condition, data = sample_meta)
  fit <- lmFit(expr, design_condition)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  do.call(rbind, lapply(names(contrast_specs), function(nm) {
    coef_name <- paste0("condition", contrast_specs[[nm]]$condition)
    tt <- topTable(fit, coef = coef_name, n = Inf, sort.by = "P")
    tt$gene_id <- rownames(tt)
    standardize_result(
      tt,
      method = "limma_trend_ComBatSeq_logCPM",
      contrast = nm,
      contrast_label = contrast_specs[[nm]]$label,
      logfc_col = "logFC",
      p_col = "P.Value",
      fdr_col = "adj.P.Val",
      stat_col = "t",
      expr_col = "AveExpr",
      model_note = "Exploratory linear-model sensitivity analysis on PCS-corrected TMM logCPM; design ~ condition."
    )
  }))
})

all_method_results <- do.call(rbind, method_results)
if (is.null(all_method_results) || nrow(all_method_results) == 0) {
  stop("No DEG sensitivity method returned results.")
}
all_method_results$method <- factor(
  all_method_results$method,
  levels = c("edgeR_QL_PCS", "edgeR_LRT_PCS", "DESeq2_PCS", "edgeR_exact_ComBatSeq",
             "DESeq2_ComBatSeq", "limma_trend_ComBatSeq_logCPM")
)
all_method_results$direction_nominal <- ifelse(all_method_results$logFC >= 0, "Up", "Down")

fold_change_grid <- c(1, 1.2, 1.5, 2)
threshold_grid <- expand.grid(
  threshold_metric = c("FDR", "PValue"),
  alpha = c(0.05, 0.10, 0.20),
  fold_change = fold_change_grid,
  stringsAsFactors = FALSE
)
threshold_grid$abs_log2fc <- log2(threshold_grid$fold_change)
threshold_grid$threshold_label <- ifelse(
  threshold_grid$fold_change <= 1,
  paste0(threshold_grid$threshold_metric, " < ", threshold_grid$alpha, ", no FC filter"),
  paste0(threshold_grid$threshold_metric, " < ", threshold_grid$alpha, ", FC >= ", threshold_grid$fold_change)
)

threshold_summary <- do.call(rbind, lapply(split(all_method_results, list(all_method_results$method, all_method_results$contrast), drop = TRUE), function(res) {
  if (nrow(res) == 0) return(data.frame())
  do.call(rbind, lapply(seq_len(nrow(threshold_grid)), function(i) {
    g <- threshold_grid[i, ]
    metric_values <- res[[g$threshold_metric]]
    sig <- !is.na(metric_values) & metric_values < g$alpha & abs(res$logFC) >= g$abs_log2fc
    data.frame(
      method = as.character(res$method[1]),
      contrast = res$contrast[1],
      contrast_label = res$contrast_label[1],
      threshold_metric = g$threshold_metric,
      alpha = g$alpha,
      fold_change = g$fold_change,
      abs_log2fc = g$abs_log2fc,
      threshold_label = g$threshold_label,
      genes_tested = nrow(res),
      up = sum(sig & res$logFC > 0, na.rm = TRUE),
      down = sum(sig & res$logFC < 0, na.rm = TRUE),
      total_deg = sum(sig, na.rm = TRUE),
      min_pvalue = suppressWarnings(min(res$PValue, na.rm = TRUE)),
      min_fdr = suppressWarnings(min(res$FDR, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  }))
}))
threshold_summary <- threshold_summary[order(threshold_summary$contrast, threshold_summary$method,
                                             threshold_summary$threshold_metric, threshold_summary$alpha,
                                             threshold_summary$fold_change), ]

key_threshold_summary <- subset(
  threshold_summary,
  (threshold_metric == "FDR" & alpha %in% c(0.05, 0.10, 0.20) & fold_change %in% c(1, 1.5, 2)) |
    (threshold_metric == "PValue" & alpha %in% c(0.05, 0.10) & fold_change %in% c(1.5, 2))
)

plot_thresholds <- subset(
  threshold_summary,
  (threshold_metric == "FDR" & alpha %in% c(0.05, 0.10, 0.20) & fold_change %in% c(1, 1.5, 2)) |
    (threshold_metric == "PValue" & alpha == 0.05 & fold_change %in% c(1.5, 2))
)
plot_thresholds$plot_label <- paste0(
  ifelse(plot_thresholds$threshold_metric == "PValue", "P", plot_thresholds$threshold_metric),
  " < ", plot_thresholds$alpha,
  "\n",
  ifelse(plot_thresholds$fold_change <= 1, "no FC filter", paste0("FC >= ", plot_thresholds$fold_change))
)
plot_thresholds$plot_label <- factor(plot_thresholds$plot_label, levels = unique(plot_thresholds$plot_label))
plot_thresholds$method <- factor(plot_thresholds$method, levels = levels(all_method_results$method))
p_heat <- ggplot(plot_thresholds, aes(x = plot_label, y = method, fill = total_deg)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = total_deg), size = 2.7, color = "grey10") +
  facet_wrap(~ contrast_label, ncol = 1) +
  scale_fill_gradient(low = "#F3F6F7", high = "#B4452D") +
  labs(
    title = "DEG calls across Shiny-style methods and thresholds",
    subtitle = "Counts are genes passing the selected p/FDR and fold-change filters; ComBat-seq methods are exploratory.",
    x = NULL,
    y = NULL,
    fill = "DEG count"
  ) +
  theme(axis.text.x = element_text(size = 8, lineheight = 0.95, hjust = 0.5)) +
  base_theme
save_gg(p_heat, "deg_method_threshold_sensitivity_heatmap", width = 12, height = 8.5)

candidate_rows <- do.call(rbind, lapply(split(all_method_results, list(all_method_results$method, all_method_results$contrast), drop = TRUE), function(res) {
  res$shiny_default <- !is.na(res$FDR) & res$FDR < 0.05 & abs(res$logFC) >= log2(1.5)
  res$publication_default <- !is.na(res$FDR) & res$FDR < 0.05 & abs(res$logFC) >= 1
  res$nominal_fc15 <- !is.na(res$PValue) & res$PValue < 0.05 & abs(res$logFC) >= log2(1.5)
  res$relaxed_fdr_fc15 <- !is.na(res$FDR) & res$FDR < 0.10 & abs(res$logFC) >= log2(1.5)
  res$selection_tier <- "top-ranked"
  res$selection_tier[res$nominal_fc15] <- "P<0.05_FC>=1.5"
  res$selection_tier[res$relaxed_fdr_fc15] <- "FDR<0.10_FC>=1.5"
  res$selection_tier[res$shiny_default] <- "FDR<0.05_FC>=1.5"
  res$selection_tier[res$publication_default] <- "FDR<0.05_FC>=2"
  res <- res[order(match(res$selection_tier, c("FDR<0.05_FC>=2", "FDR<0.05_FC>=1.5", "FDR<0.10_FC>=1.5", "P<0.05_FC>=1.5", "top-ranked")),
                   res$FDR, res$PValue, -abs(res$logFC), na.last = TRUE), ]
  head(res, 25)
}))
candidate_cols <- intersect(
  c("method", "contrast_label", "selection_tier", "gene_id", "gene_label", "SYMBOL", "GENENAME",
    "logFC", "statistic", "PValue", "FDR", "average_expression", "model_note"),
  names(candidate_rows)
)
top_multimethod_candidates <- candidate_rows[, candidate_cols, drop = FALSE]

primary_diagnostics <- do.call(rbind, lapply(names(s2$edgeR_results), function(nm) {
  res <- s2$edgeR_results[[nm]]
  top <- res[order(res$PValue, res$FDR, -abs(res$logFC), na.last = TRUE), ][1, ]
  data.frame(
    contrast = nm,
    contrast_label = contrast_specs[[nm]]$label,
    genes_tested = nrow(res),
    min_pvalue = min(res$PValue, na.rm = TRUE),
    min_fdr = min(res$FDR, na.rm = TRUE),
    raw_p_lt_0_05 = sum(res$PValue < 0.05, na.rm = TRUE),
    raw_p_lt_0_10 = sum(res$PValue < 0.10, na.rm = TRUE),
    fdr_lt_0_05 = sum(res$FDR < 0.05, na.rm = TRUE),
    abs_log2fc_ge_log2_1_5 = sum(abs(res$logFC) >= log2(1.5), na.rm = TRUE),
    abs_log2fc_ge_1 = sum(abs(res$logFC) >= 1, na.rm = TRUE),
    raw_p_lt_0_05_and_fc_ge_1_5 = sum(res$PValue < 0.05 & abs(res$logFC) >= log2(1.5), na.rm = TRUE),
    raw_p_lt_0_05_and_fc_ge_2 = sum(res$PValue < 0.05 & abs(res$logFC) >= 1, na.rm = TRUE),
    top_gene = top$gene_label,
    top_gene_logFC = top$logFC,
    top_gene_PValue = top$PValue,
    top_gene_FDR = top$FDR,
    stringsAsFactors = FALSE
  )
}))

bcv <- sqrt(s2$y$tagwise.dispersion)
dispersion_summary <- data.frame(
  metric = c("median_BCV", "IQR_BCV", "90th_percentile_BCV", "median_edgeR_residual_df", "edgeR_prior_df"),
  value = c(
    median(bcv, na.rm = TRUE),
    IQR(bcv, na.rm = TRUE),
    as.numeric(quantile(bcv, 0.90, na.rm = TRUE)),
    median(s2$fit$df.residual, na.rm = TRUE),
    median(s2$fit$df.prior, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

gene_vars <- apply(s1$raw_logcpm, 1, var)
vp_genes <- names(sort(gene_vars, decreasing = TRUE))[seq_len(min(5000, length(gene_vars)))]
eta_mat <- t(vapply(vp_genes, function(g) {
  dat <- data.frame(expr = as.numeric(s1$raw_logcpm[g, sample_meta$sample_id]), pcs = sample_meta$pcs, condition = sample_meta$condition)
  a <- anova(lm(expr ~ pcs + condition, data = dat))
  ss <- setNames(a[, "Sum Sq"], rownames(a))
  total <- sum(ss, na.rm = TRUE)
  c(
    eta2_pcs = unname(ss["pcs"] / total),
    eta2_condition = unname(ss["condition"] / total),
    eta2_residual = unname(ss["Residuals"] / total)
  )
}, numeric(3)))
variance_partition_summary <- data.frame(
  metric = c("top_variable_genes_evaluated", "median_eta2_PCS", "median_eta2_condition", "median_eta2_residual",
             "genes_with_PCS_eta2_gt_condition_eta2"),
  value = c(
    nrow(eta_mat),
    median(eta_mat[, "eta2_pcs"], na.rm = TRUE),
    median(eta_mat[, "eta2_condition"], na.rm = TRUE),
    median(eta_mat[, "eta2_residual"], na.rm = TRUE),
    sum(eta_mat[, "eta2_pcs"] > eta_mat[, "eta2_condition"], na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

best_relaxed <- threshold_summary[order(-threshold_summary$total_deg, threshold_summary$threshold_metric,
                                        threshold_summary$alpha, -threshold_summary$fold_change), ][1, ]
reason_table <- data.frame(
  reason = c(
    "Limited paired sample degrees of freedom",
    "PCS block heterogeneity is large relative to treatment signal",
    "Multiple-testing correction is stringent with more than twenty thousand genes",
    "Observed single-gene effects are modest under the primary count model",
    "Publication and Shiny-style FDR thresholds are more conservative than nominal P-value screens",
    "Exploratory batch-corrected count tests are useful for candidate finding, not primary inference",
    "Unmodeled upstream technical variation cannot be ruled out"
  ),
  evidence = c(
    paste0(length(unique(sample_meta$pcs)), " PCS blocks with one sample per condition per block; median edgeR residual df = ",
           round(dispersion_summary$value[dispersion_summary$metric == "median_edgeR_residual_df"], 2), "."),
    paste0("Across the top variable genes, median eta2 for PCS = ",
           round(variance_partition_summary$value[variance_partition_summary$metric == "median_eta2_PCS"], 3),
           " versus median eta2 for condition = ",
           round(variance_partition_summary$value[variance_partition_summary$metric == "median_eta2_condition"], 3), "."),
    paste0(format(primary_diagnostics$genes_tested[1], big.mark = ","), " genes are tested per contrast; primary minimum FDR values are ",
           paste(paste0(primary_diagnostics$contrast, ": ", signif(primary_diagnostics$min_fdr, 3)), collapse = "; "), "."),
    paste0("Primary edgeR raw P<0.05 and FC>=1.5 counts are ",
           paste(paste0(primary_diagnostics$contrast, ": ", primary_diagnostics$raw_p_lt_0_05_and_fc_ge_1_5), collapse = "; "),
           "."),
    paste0("The largest DEG count in the sensitivity grid is ", best_relaxed$total_deg, " for ",
           best_relaxed$method, ", ", best_relaxed$contrast, " at ", best_relaxed$threshold_label, "."),
    "ComBat-seq condition-preserving counts are used downstream to mirror the Shiny app's batch-adjusted exact/DESeq2 style checks, but the primary DEG table remains the raw-count model with PCS in the design.",
    "No FASTQ/BAM-level QC, RNA integrity, strandedness, alignment rate, or library-prep covariates were supplied to test additional sources of variation."
  ),
  interpretation = c(
    "Small n limits power after fitting five PCS block terms and two treatment coefficients.",
    "A strong paired donor/block structure can absorb variance that might otherwise appear as treatment signal.",
    "Many nominally small P values can become non-significant after Benjamini-Hochberg correction.",
    "Single genes may trend consistently but remain below the effect-size and uncertainty cutoffs.",
    "Relaxed or nominal thresholds can suggest candidates, but they should be presented as exploratory.",
    "Candidate genes shared across primary and exploratory methods are more credible than genes unique to a corrected-count screen.",
    "Additional covariates or upstream QC failures could mask treatment effects or inflate dispersion."
  ),
  stringsAsFactors = FALSE
)

write.csv(all_method_results, file.path(table_dir, "DEG_method_sensitivity_all_gene_results.csv"), row.names = FALSE)
write.csv(threshold_summary, file.path(table_dir, "DEG_method_threshold_sensitivity_summary.csv"), row.names = FALSE)
write.csv(key_threshold_summary, file.path(table_dir, "DEG_method_threshold_key_summary.csv"), row.names = FALSE)
write.csv(top_multimethod_candidates, file.path(table_dir, "DEG_multimethod_top_candidate_genes.csv"), row.names = FALSE)
write.csv(primary_diagnostics, file.path(table_dir, "DEG_primary_edgeR_no_deg_diagnostics.csv"), row.names = FALSE)
write.csv(dispersion_summary, file.path(table_dir, "DEG_edgeR_dispersion_power_diagnostics.csv"), row.names = FALSE)
write.csv(variance_partition_summary, file.path(table_dir, "DEG_variance_partition_top_variable_genes.csv"), row.names = FALSE)
write.csv(reason_table, file.path(table_dir, "possible_reasons_for_few_or_no_DEGs.csv"), row.names = FALSE)
write.csv(method_status, file.path(table_dir, "DEG_method_sensitivity_status.csv"), row.names = FALSE)

reason_md <- c(
  "# Possible reasons for few or no high-confidence DEGs",
  "",
  "This note is generated from the cached multi-method sensitivity stage. It should be read together with the primary edgeR tables and the method/threshold summary grid.",
  "",
  apply(reason_table, 1, function(x) {
    paste0("## ", x[["reason"]], "\n\nEvidence: ", x[["evidence"]], "\n\nInterpretation: ", x[["interpretation"]], "\n")
  })
)
writeLines(reason_md, file.path(report_dir, "possible_reasons_for_few_or_no_DEGs.md"))

stage_data <- list(
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  sensitivity_model = "Shiny-app-style DEG sensitivity: edgeR QL/LRT with PCS covariate, edgeR exact on ComBat-seq counts, DESeq2 with PCS, DESeq2 on ComBat-seq counts, limma-trend on PCS-corrected logCPM.",
  all_method_results = all_method_results,
  threshold_summary = threshold_summary,
  key_threshold_summary = key_threshold_summary,
  top_multimethod_candidates = top_multimethod_candidates,
  primary_diagnostics = primary_diagnostics,
  dispersion_summary = dispersion_summary,
  variance_partition_summary = variance_partition_summary,
  possible_reasons = reason_table,
  method_status = method_status
)
saveRDS(stage_data, stage_rds)
write_stage_stamp(stage_name)
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_05_method_threshold_sensitivity.txt"))
message("Stage 5 complete: ", stage_rds)
