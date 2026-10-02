# Scientific workflow derived from Revision/RNAseq PCS met inhibitors/results/scripts/02_deg_edger.R
# Run from the repository root with METABOLIC_INPUT_DIR set to the external project data.
source("data_preprocessing/common/config.R")
script_dir <- {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)) else getwd()
}
source(file.path(script_dir, "00_config.R"))

stage_name <- "02_deg_edger"
stage_rds <- file.path(cache_dir, paste0(stage_name, ".rds"))
stage1_rds <- file.path(cache_dir, "01_qc_batch_correction.rds")
deps <- c(stage1_rds, file.path(script_dir, "00_config.R"), file.path(script_dir, "02_deg_edger.R"))

if (stage_is_fresh(stage_rds, deps)) {
  message("Using cached ", stage_name, ": ", stage_rds)
  quit(save = "no", status = 0)
}

message("Stage 2: edgeR quasi-likelihood DEG analysis with PCS as batch covariate")
s1 <- read_stage(stage1_rds)
y <- s1$y
design <- s1$design_batch_condition
annotation <- s1$filtered_annotation
sample_meta <- s1$sample_meta

y <- estimateDisp(y, design = design, robust = TRUE)
fit <- glmQLFit(y, design = design, robust = TRUE)

condition_means <- data.frame(
  gene_id = rownames(y$counts),
  mean_count_C = rowMeans(y$counts[, sample_meta$condition == "C", drop = FALSE]),
  mean_count_2DG = rowMeans(y$counts[, sample_meta$condition == "2DG", drop = FALSE]),
  mean_count_LDHA = rowMeans(y$counts[, sample_meta$condition == "LDHA", drop = FALSE]),
  mean_batch_corrected_count_C = rowMeans(s1$batch_corrected_pseudo_counts[, sample_meta$condition == "C", drop = FALSE]),
  mean_batch_corrected_count_2DG = rowMeans(s1$batch_corrected_pseudo_counts[, sample_meta$condition == "2DG", drop = FALSE]),
  mean_batch_corrected_count_LDHA = rowMeans(s1$batch_corrected_pseudo_counts[, sample_meta$condition == "LDHA", drop = FALSE]),
  stringsAsFactors = FALSE
)

edgeR_results <- list()
for (nm in names(contrast_specs)) {
  spec <- contrast_specs[[nm]]
  qlf <- glmQLFTest(fit, coef = spec$coef)
  tt <- as.data.frame(topTags(qlf, n = Inf, sort.by = "PValue"))
  tt$gene_id <- rownames(tt)
  tt <- merge(tt, annotation, by = "gene_id", all.x = TRUE, sort = FALSE)
  tt <- merge(tt, condition_means, by = "gene_id", all.x = TRUE, sort = FALSE)
  tt <- tt[match(rownames(topTags(qlf, n = Inf, sort.by = "PValue")$table), tt$gene_id), ]
  tt$contrast <- nm
  tt$contrast_label <- spec$label
  tt$direction <- classify_deg(tt)
  tt$neg_log10_fdr <- -log10(pmax(tt$FDR, .Machine$double.xmin))
  tt$rank_metric_signed_sqrt_F <- sign(tt$logFC) * sqrt(pmax(tt$F, 0))
  edgeR_results[[nm]] <- tt
  write.csv(tt, file.path(table_dir, paste0("edgeR_QLF_PCS_batch_", nm, "_all_genes.csv")), row.names = FALSE)
  write.csv(tt[tt$FDR < 0.05 & abs(tt$logFC) >= 1, ],
            file.path(table_dir, paste0("edgeR_QLF_PCS_batch_", nm, "_DEG_FDR05_abslog2FC1.csv")), row.names = FALSE)
}

deg_summary <- do.call(rbind, lapply(names(edgeR_results), function(nm) {
  res <- edgeR_results[[nm]]
  data.frame(
    contrast = nm,
    contrast_label = contrast_specs[[nm]]$label,
    genes_tested = nrow(res),
    fdr_0_05 = sum(res$FDR < 0.05, na.rm = TRUE),
    up_fdr_0_05_abslog2fc_1 = sum(res$FDR < 0.05 & res$logFC >= 1, na.rm = TRUE),
    down_fdr_0_05_abslog2fc_1 = sum(res$FDR < 0.05 & res$logFC <= -1, na.rm = TRUE),
    deg_fdr_0_05_abslog2fc_1 = sum(res$FDR < 0.05 & abs(res$logFC) >= 1, na.rm = TRUE)
  )
}))
write.csv(deg_summary, file.path(table_dir, "edgeR_QLF_PCS_batch_contrast_summary.csv"), row.names = FALSE)

plot_volcano <- function(res, stem, title) {
  labels <- res[res$FDR < 0.05 & abs(res$logFC) >= 1, ]
  labels <- labels[order(labels$FDR, -abs(labels$logFC)), ]
  labels <- head(labels, 16)
  p <- ggplot(res, aes(x = logFC, y = neg_log10_fdr, color = direction)) +
    geom_point(size = 1.1, alpha = 0.75, stroke = 0) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey45", linewidth = 0.35) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey45", linewidth = 0.35) +
    ggrepel::geom_text_repel(
      data = labels,
      aes(label = gene_label),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.35,
      min.segment.length = 0,
      show.legend = FALSE
    ) +
    scale_color_manual(values = deg_colors, drop = TRUE) +
    labs(
      title = title,
      subtitle = "edgeR QL GLM with PCS batch covariate; dashed lines mark FDR 0.05 and |log2FC| 1",
      x = "log2 fold-change",
      y = "-log10 adjusted P value",
      color = "Class"
    ) +
    base_theme
  save_gg(p, stem, width = 7, height = 5.8)
}

plot_ma <- function(res, stem, title) {
  p <- ggplot(res, aes(x = logCPM, y = logFC, color = direction)) +
    geom_hline(yintercept = 0, color = "grey45", linewidth = 0.35) +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "grey65", linewidth = 0.3) +
    geom_point(size = 1.0, alpha = 0.72, stroke = 0) +
    scale_color_manual(values = deg_colors, drop = TRUE) +
    labs(
      title = title,
      subtitle = "Average expression is edgeR logCPM; dashed lines mark |log2FC| 1",
      x = "Average logCPM",
      y = "log2 fold-change",
      color = "Class"
    ) +
    base_theme
  save_gg(p, stem, width = 7, height = 5.2)
}

for (nm in names(edgeR_results)) {
  plot_volcano(edgeR_results[[nm]], paste0("deg_edgeR_volcano_", nm), paste0("Volcano plot: ", contrast_specs[[nm]]$label))
  plot_ma(edgeR_results[[nm]], paste0("deg_edgeR_ma_", nm), paste0("MA plot: ", contrast_specs[[nm]]$label))
}

deg_long <- do.call(rbind, lapply(names(edgeR_results), function(nm) {
  res <- edgeR_results[[nm]]
  data.frame(contrast = contrast_specs[[nm]]$label, direction = as.character(res$direction), stringsAsFactors = FALSE)
}))
deg_bar_df <- as.data.frame(table(deg_long$contrast, deg_long$direction), stringsAsFactors = FALSE)
colnames(deg_bar_df) <- c("contrast", "direction", "n")
deg_bar_df <- deg_bar_df[deg_bar_df$direction %in% c("Up", "Down", "FDR only"), ]
if (nrow(deg_bar_df) == 0) {
  deg_bar_df <- data.frame(contrast = vapply(contrast_specs, `[[`, character(1), "label"), direction = "No DEG", n = 0)
}
p_deg_bar <- ggplot(deg_bar_df, aes(x = contrast, y = n, fill = direction)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65, color = "grey25", linewidth = 0.2) +
  geom_text(aes(label = n), position = position_dodge(width = 0.75), vjust = -0.35, size = 3) +
  scale_fill_manual(values = c(deg_colors, `No DEG` = "#BDBDBD"), drop = TRUE) +
  labs(
    title = "edgeR differential expression burden by treatment",
    subtitle = "Up/down require FDR < 0.05 and |log2FC| >= 1",
    x = NULL,
    y = "Number of genes",
    fill = "Class"
  ) +
  coord_cartesian(ylim = c(0, max(deg_bar_df$n, 1) * 1.18)) +
  base_theme
save_gg(p_deg_bar, "deg_edgeR_summary_bar", width = 7, height = 4.8)

top_gene_ids <- unique(unlist(lapply(edgeR_results, function(res) {
  sig <- res[res$FDR < 0.05 & abs(res$logFC) >= 1, ]
  if (nrow(sig) < 20) sig <- res
  head(sig$gene_id[order(sig$FDR, -abs(sig$logFC))], 50)
})))
top_gene_ids <- intersect(top_gene_ids, rownames(s1$batch_corrected_logcpm))
heat_mat <- s1$batch_corrected_logcpm[top_gene_ids, sample_meta$sample_id, drop = FALSE]
heat_labels <- annotation$gene_label[match(rownames(heat_mat), annotation$gene_id)]
heat_labels <- make.unique(ifelse(is.na(heat_labels) | heat_labels == "", rownames(heat_mat), heat_labels))
rownames(heat_mat) <- heat_labels
ann_col_heat <- data.frame(
  PCS = sample_meta$pcs,
  Condition = sample_meta$condition
)
rownames(ann_col_heat) <- sample_meta$sample_id
ann_colors <- list(
  Condition = condition_colors,
  PCS = setNames(brewer.pal(max(3, length(levels(sample_meta$pcs))), "Set2")[seq_along(levels(sample_meta$pcs))],
                 levels(sample_meta$pcs))
)
pheatmap(
  heat_mat,
  scale = "row",
  annotation_col = ann_col_heat,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("#2166ac", "#f7f7f7", "#b35806"))(100),
  border_color = NA,
  fontsize_row = ifelse(nrow(heat_mat) > 80, 5, 7),
  fontsize_col = 8,
  main = "Top edgeR genes on PCS-corrected logCPM",
  filename = file.path(figure_dir, "deg_edgeR_top_genes_PCS_corrected_heatmap.png"),
  width = 8.5,
  height = max(6.5, min(14, 2 + nrow(heat_mat) * 0.12))
)
pheatmap(
  heat_mat,
  scale = "row",
  annotation_col = ann_col_heat,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("#2166ac", "#f7f7f7", "#b35806"))(100),
  border_color = NA,
  fontsize_row = ifelse(nrow(heat_mat) > 80, 5, 7),
  fontsize_col = 8,
  main = "Top edgeR genes on PCS-corrected logCPM",
  filename = file.path(figure_dir, "deg_edgeR_top_genes_PCS_corrected_heatmap.pdf"),
  width = 8.5,
  height = max(6.5, min(14, 2 + nrow(heat_mat) * 0.12))
)

top_genes_summary <- do.call(rbind, lapply(names(edgeR_results), function(nm) {
  res <- edgeR_results[[nm]]
  sig <- res[res$FDR < 0.05 & abs(res$logFC) >= 1, ]
  if (nrow(sig) == 0) sig <- head(res[order(res$FDR), ], 10)
  sig <- sig[order(sig$FDR, -abs(sig$logFC)), ]
  head(sig[, c("contrast_label", "gene_id", "gene_label", "SYMBOL", "GENENAME", "logFC", "logCPM", "F", "PValue", "FDR", "direction")], 12)
}))
write.csv(top_genes_summary, file.path(table_dir, "edgeR_top_genes_summary_for_report.csv"), row.names = FALSE)

stage_data <- list(
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  edgeR_model = "edgeR glmQLFit/glmQLFTest with TMM normalization and design ~ PCS + condition",
  y = y,
  design = design,
  fit = fit,
  edgeR_results = edgeR_results,
  deg_summary = deg_summary,
  top_genes_summary = top_genes_summary
)
saveRDS(stage_data, stage_rds)
write_stage_stamp(stage_name)
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_02_deg_edger.txt"))
message("Stage 2 complete: ", stage_rds)
