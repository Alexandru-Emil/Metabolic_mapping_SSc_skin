# Scientific workflow derived from Revision/RNAseq PCS met inhibitors/results/scripts/01_qc_batch_correction.R
# Run from the repository root with METABOLIC_INPUT_DIR set to the external project data.
source("data_preprocessing/common/config.R")
script_dir <- {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)) else getwd()
}
source(file.path(script_dir, "00_config.R"))

stage_name <- "01_qc_batch_correction"
stage_rds <- file.path(cache_dir, paste0(stage_name, ".rds"))
deps <- c(input_file, file.path(script_dir, "00_config.R"), file.path(script_dir, "01_qc_batch_correction.R"))

if (stage_is_fresh(stage_rds, deps)) {
  message("Using cached ", stage_name, ": ", stage_rds)
  quit(save = "no", status = 0)
}

message("Stage 1: QC, technical-lane aggregation, and PCS batch correction")
if (!file.exists(input_file)) stop("Input file not found: ", input_file)

first_line <- readLines(input_file, n = 1, warn = FALSE)
lane_names <- strsplit(first_line, "\t", fixed = TRUE)[[1]]
raw_counts_df <- read.delim(
  input_file,
  header = FALSE,
  skip = 1,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (ncol(raw_counts_df) != length(lane_names) + 1) {
  stop("Unexpected table shape. Header has ", length(lane_names),
       " lane names, data rows have ", ncol(raw_counts_df), " columns.")
}

colnames(raw_counts_df) <- c("gene_id", lane_names)
gene_ids <- raw_counts_df$gene_id
lane_counts <- as.matrix(raw_counts_df[, lane_names, drop = FALSE])
storage.mode(lane_counts) <- "numeric"
rownames(lane_counts) <- gene_ids

if (anyDuplicated(rownames(lane_counts))) {
  lane_counts <- rowsum(lane_counts, group = rownames(lane_counts), reorder = FALSE)
  gene_ids <- rownames(lane_counts)
}

lane_meta <- do.call(rbind, lapply(lane_names, parse_lane_name))
lane_meta$sample_id <- paste(lane_meta$pcs, lane_meta$condition, sep = "_")
lane_meta$lane_library_size <- colSums(lane_counts[, lane_meta$lane_name, drop = FALSE])
lane_meta$condition <- factor(lane_meta$condition, levels = c("C", "2DG", "LDHA"))
lane_meta$condition_label <- condition_labels[as.character(lane_meta$condition)]

sample_meta <- unique(lane_meta[, c("sample_id", "pcs", "condition", "condition_label")])
sample_meta$pcs_number <- as.integer(sub("^PCS", "", sample_meta$pcs))
sample_meta <- sample_meta[order(sample_meta$pcs_number, match(sample_meta$condition, c("C", "2DG", "LDHA"))), ]
rownames(sample_meta) <- sample_meta$sample_id
sample_meta$pcs <- factor(sample_meta$pcs, levels = paste0("PCS", sort(unique(sample_meta$pcs_number))))
sample_meta$condition <- factor(sample_meta$condition, levels = c("C", "2DG", "LDHA"))
lane_count_by_sample <- table(lane_meta$sample_id)
sample_meta$n_technical_lanes <- as.integer(lane_count_by_sample[sample_meta$sample_id])

aggregated_counts <- sapply(sample_meta$sample_id, function(sid) {
  rowSums(lane_counts[, lane_meta$sample_id == sid, drop = FALSE])
})
aggregated_counts <- as.matrix(aggregated_counts)
storage.mode(aggregated_counts) <- "numeric"
rownames(aggregated_counts) <- gene_ids
sample_meta$library_size <- colSums(aggregated_counts[, sample_meta$sample_id, drop = FALSE])

rounded_counts <- round(aggregated_counts)
rounding_delta <- aggregated_counts - rounded_counts
rounding_summary <- data.frame(
  metric = c("max_abs_rounding_delta", "mean_abs_rounding_delta", "fraction_values_changed_by_rounding"),
  value = c(max(abs(rounding_delta)), mean(abs(rounding_delta)), mean(abs(rounding_delta) > 0))
)

annotation <- annotate_genes(rownames(aggregated_counts))
sample_meta_model <- sample_meta
sample_meta_model$condition <- relevel(sample_meta_model$condition, ref = "C")
design_batch_condition <- model.matrix(~ pcs + condition, data = sample_meta_model)
rownames(design_batch_condition) <- sample_meta_model$sample_id
design_condition <- model.matrix(~ condition, data = sample_meta_model)
rownames(design_condition) <- sample_meta_model$sample_id

y0 <- DGEList(counts = rounded_counts[, sample_meta_model$sample_id, drop = FALSE])
nonzero <- rowSums(y0$counts) > 0
y0 <- y0[nonzero, , keep.lib.sizes = FALSE]
keep <- filterByExpr(y0, design = design_batch_condition)
y <- y0[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")

filtered_counts <- y$counts
filtered_annotation <- annotation[match(rownames(filtered_counts), annotation$gene_id), ]
raw_logcpm <- cpm(y, log = TRUE, prior.count = 2)

message("Running ComBat-seq with PCS as batch and condition preserved.")
batch_corrected_counts <- tryCatch({
  bc <- ComBat_seq(
    counts = as.matrix(filtered_counts),
    batch = sample_meta_model$pcs,
    group = sample_meta_model$condition
  )
  storage.mode(bc) <- "numeric"
  bc
}, error = function(e) {
  writeLines(conditionMessage(e), file.path(log_dir, "ComBat_seq_error.txt"))
  message("ComBat-seq failed; falling back to removeBatchEffect on logCPM.")
  NULL
})

if (!is.null(batch_corrected_counts)) {
  y_bc <- DGEList(counts = batch_corrected_counts)
  y_bc <- calcNormFactors(y_bc, method = "TMM")
  batch_corrected_logcpm <- cpm(y_bc, log = TRUE, prior.count = 2)
  batch_correction_method <- "sva::ComBat_seq on rounded, filtered, lane-summed counts; condition preserved"
  batch_corrected_pseudo_counts <- batch_corrected_counts
} else {
  batch_corrected_logcpm <- removeBatchEffect(raw_logcpm, batch = sample_meta_model$pcs, design = design_condition)
  batch_correction_method <- "limma::removeBatchEffect on TMM logCPM fallback; condition preserved"
  mean_lib <- mean(sample_meta_model$library_size)
  batch_corrected_pseudo_counts <- 2^batch_corrected_logcpm * mean_lib / 1e6
}

filtering_summary <- data.frame(
  metric = c("input_genes", "nonzero_genes", "genes_after_filterByExpr", "biological_samples", "sequencing_lanes"),
  value = c(nrow(aggregated_counts), sum(nonzero), nrow(filtered_counts), ncol(aggregated_counts), ncol(lane_counts))
)

write.csv(lane_meta, file.path(table_dir, "lane_metadata.csv"), row.names = FALSE)
write.csv(sample_meta, file.path(table_dir, "sample_metadata.csv"), row.names = FALSE)
write.csv(annotation, file.path(table_dir, "gene_annotation.csv"), row.names = FALSE)
write.csv(aggregated_counts, file.path(table_dir, "technical_lane_summed_counts_fractional.csv"), row.names = TRUE)
write.csv(rounded_counts, file.path(table_dir, "technical_lane_summed_counts_rounded_for_edgeR.csv"), row.names = TRUE)
write.csv(rounding_summary, file.path(table_dir, "rounding_summary_for_edgeR.csv"), row.names = FALSE)
write.csv(filtering_summary, file.path(table_dir, "filtering_summary_edgeR.csv"), row.names = FALSE)
write.csv(raw_logcpm, file.path(table_dir, "tmm_logCPM_before_PCS_batch_correction.csv"), row.names = TRUE)
write.csv(batch_corrected_logcpm, file.path(table_dir, "PCS_batch_corrected_logCPM.csv"), row.names = TRUE)
write.csv(batch_corrected_pseudo_counts, file.path(table_dir, "PCS_batch_corrected_counts_for_downstream_visualization.csv"), row.names = TRUE)

condition_colors_plot <- condition_colors
names(condition_colors_plot) <- c("C", "2DG", "LDHA")
lib_df <- data.frame(
  sample_id = sample_meta_model$sample_id,
  pcs = sample_meta_model$pcs,
  condition = sample_meta_model$condition,
  condition_label = sample_meta_model$condition_label,
  n_technical_lanes = sample_meta_model$n_technical_lanes,
  library_size_millions = sample_meta_model$library_size / 1e6
)
p_lib <- ggplot(lib_df, aes(x = sample_id, y = library_size_millions, fill = condition)) +
  geom_col(width = 0.72, color = "grey25", linewidth = 0.2) +
  geom_text(aes(label = n_technical_lanes), vjust = -0.4, size = 3, color = "grey20") +
  scale_fill_manual(values = condition_colors_plot, labels = condition_labels) +
  labs(
    title = "Library sizes after technical lane consolidation",
    subtitle = "Numbers above bars show technical lanes summed per biological sample",
    x = NULL,
    y = "Summed count depth (millions)",
    fill = "Condition"
  ) +
  coord_cartesian(ylim = c(0, max(lib_df$library_size_millions) * 1.14)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  base_theme
save_gg(p_lib, "qc_library_sizes_edgeR_inputs", width = 8, height = 4.8)

make_pca_df <- function(mat, label) {
  gene_vars <- apply(mat, 1, var)
  pca_genes <- names(sort(gene_vars, decreasing = TRUE))[seq_len(min(5000, length(gene_vars)))]
  pca <- prcomp(t(mat[pca_genes, , drop = FALSE]), center = TRUE, scale. = FALSE)
  pve <- (pca$sdev^2) / sum(pca$sdev^2)
  out <- cbind(sample_meta_model, as.data.frame(pca$x[, 1:2, drop = FALSE]))
  out$matrix <- label
  out$xlab <- paste0("PC1 (", round(100 * pve[1], 1), "%)")
  out$ylab <- paste0("PC2 (", round(100 * pve[2], 1), "%)")
  attr(out, "pve") <- pve
  out
}

pca_raw <- make_pca_df(raw_logcpm, "Before PCS correction")
pca_bc <- make_pca_df(batch_corrected_logcpm, "After PCS correction")
pca_df <- rbind(pca_raw, pca_bc)
pca_df$matrix <- factor(pca_df$matrix, levels = c("Before PCS correction", "After PCS correction"))
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, shape = pcs, label = sample_id)) +
  geom_hline(yintercept = 0, color = "grey88", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey88", linewidth = 0.3) +
  geom_point(size = 3.1, stroke = 0.7) +
  ggrepel::geom_text_repel(size = 2.7, max.overlaps = Inf, box.padding = 0.3,
                           min.segment.length = 0, show.legend = FALSE) +
  facet_wrap(~ matrix, scales = "free") +
  scale_color_manual(values = condition_colors_plot, labels = condition_labels) +
  labs(
    title = "PCA before and after PCS batch correction",
    subtitle = "Top 5,000 most variable retained genes; ComBat-seq preserves treatment groups",
    x = "PC1",
    y = "PC2",
    color = "Condition",
    shape = "PCS block"
  ) +
  base_theme
save_gg(p_pca, "qc_pca_before_after_PCS_batch_correction", width = 10, height = 5.7)

sample_cor <- cor(batch_corrected_logcpm, method = "pearson")
ann_col <- data.frame(
  PCS = sample_meta_model$pcs,
  Condition = sample_meta_model$condition
)
rownames(ann_col) <- sample_meta_model$sample_id
ann_colors <- list(
  Condition = condition_colors_plot,
  PCS = setNames(brewer.pal(max(3, length(levels(sample_meta_model$pcs))), "Set2")[seq_along(levels(sample_meta_model$pcs))],
                 levels(sample_meta_model$pcs))
)
pheatmap(
  sample_cor,
  annotation_col = ann_col,
  annotation_row = ann_col,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("#f7fbff", "#6baed6", "#08306b"))(100),
  border_color = NA,
  fontsize = 8,
  main = "Sample correlation after PCS correction",
  filename = file.path(figure_dir, "qc_sample_correlation_after_PCS_batch_correction.png"),
  width = 7.5,
  height = 6.8
)
pheatmap(
  sample_cor,
  annotation_col = ann_col,
  annotation_row = ann_col,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("#f7fbff", "#6baed6", "#08306b"))(100),
  border_color = NA,
  fontsize = 8,
  main = "Sample correlation after PCS correction",
  filename = file.path(figure_dir, "qc_sample_correlation_after_PCS_batch_correction.pdf"),
  width = 7.5,
  height = 6.8
)

stage_data <- list(
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  input_file = input_file,
  batch_correction_method = batch_correction_method,
  lane_meta = lane_meta,
  sample_meta = sample_meta_model,
  annotation = annotation,
  filtered_annotation = filtered_annotation,
  aggregated_counts_fractional = aggregated_counts,
  rounded_counts = rounded_counts,
  filtered_counts = filtered_counts,
  raw_logcpm = raw_logcpm,
  batch_corrected_logcpm = batch_corrected_logcpm,
  batch_corrected_pseudo_counts = batch_corrected_pseudo_counts,
  y = y,
  design_batch_condition = design_batch_condition,
  filtering_summary = filtering_summary,
  rounding_summary = rounding_summary
)
saveRDS(stage_data, stage_rds)
write_stage_stamp(stage_name)
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_01_qc_batch_correction.txt"))
message("Stage 1 complete: ", stage_rds)
