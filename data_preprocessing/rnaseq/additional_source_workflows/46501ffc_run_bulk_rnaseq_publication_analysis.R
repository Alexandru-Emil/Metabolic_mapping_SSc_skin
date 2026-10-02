source("data_preprocessing/common/config.R")
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(openxlsx)
  library(rmarkdown)
  library(knitr)
})

set.seed(42)

root_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
input_file <- file.path(root_dir, "inputs", "Galaxy452-gene_level_values.tabular")
analysis_dir <- file.path(root_dir, "bulk_rnaseq_2dg_ldha_analysis")
table_dir <- file.path(analysis_dir, "tables")
figure_dir <- file.path(analysis_dir, "figures")
object_dir <- file.path(analysis_dir, "objects")
report_dir <- file.path(analysis_dir, "report")
log_dir <- file.path(analysis_dir, "logs")

for (d in c(analysis_dir, table_dir, figure_dir, object_dir, report_dir, log_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

message("Reading count table: ", input_file)
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

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
  message("Duplicate gene IDs detected; summing duplicated rows.")
  lane_counts <- rowsum(lane_counts, group = rownames(lane_counts), reorder = FALSE)
  gene_ids <- rownames(lane_counts)
}

parse_lane <- function(x) {
  m <- regexec("^(PCS[0-9]+)_(C|2DG|LDHA)_(.+)$", x)
  hit <- regmatches(x, m)[[1]]
  if (length(hit) != 4) {
    stop("Could not parse sample name: ", x)
  }
  data.frame(
    lane_name = x,
    pcs = hit[2],
    condition = hit[3],
    technical_label = hit[4],
    stringsAsFactors = FALSE
  )
}

lane_meta <- do.call(rbind, lapply(lane_names, parse_lane))
lane_meta$sample_id <- paste(lane_meta$pcs, lane_meta$condition, sep = "_")
lane_meta$lane_library_size <- colSums(lane_counts[, lane_meta$lane_name, drop = FALSE])
lane_meta$condition <- factor(lane_meta$condition, levels = c("C", "2DG", "LDHA"))
lane_meta$condition_label <- c(C = "Control", `2DG` = "2DG + OX", LDHA = "LDHA + OX")[as.character(lane_meta$condition)]

sample_meta <- unique(lane_meta[, c("sample_id", "pcs", "condition", "condition_label")])
sample_meta$pcs_number <- as.integer(sub("^PCS", "", sample_meta$pcs))
sample_meta <- sample_meta[order(sample_meta$pcs_number, match(sample_meta$condition, c("C", "2DG", "LDHA"))), ]
rownames(sample_meta) <- sample_meta$sample_id
sample_meta$pcs <- factor(sample_meta$pcs, levels = paste0("PCS", sort(unique(sample_meta$pcs_number))))
sample_meta$condition <- factor(sample_meta$condition, levels = c("C", "2DG", "LDHA"))

lane_count_by_sample <- table(lane_meta$sample_id)
sample_meta$n_technical_lanes <- as.integer(lane_count_by_sample[sample_meta$sample_id])

agg_counts <- sapply(sample_meta$sample_id, function(sid) {
  rowSums(lane_counts[, lane_meta$sample_id == sid, drop = FALSE])
})
agg_counts <- as.matrix(agg_counts)
storage.mode(agg_counts) <- "numeric"
rownames(agg_counts) <- gene_ids
sample_meta$library_size <- colSums(agg_counts[, sample_meta$sample_id, drop = FALSE])

write.csv(lane_meta, file.path(table_dir, "lane_metadata.csv"), row.names = FALSE)
write.csv(sample_meta, file.path(table_dir, "sample_metadata.csv"), row.names = FALSE)
write.csv(agg_counts, file.path(table_dir, "technical_lane_summed_counts.csv"), row.names = TRUE)

message("Annotating Ensembl IDs.")
ensembl_base <- sub("\\..*$", "", rownames(agg_counts))
annotation <- data.frame(
  gene_id = rownames(agg_counts),
  ensembl_id = ensembl_base,
  stringsAsFactors = FALSE
)

ann_map <- suppressMessages(
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(ensembl_base),
    keytype = "ENSEMBL",
    columns = c("SYMBOL", "GENENAME", "ENTREZID")
  )
)
ann_map <- ann_map[order(is.na(ann_map$SYMBOL), is.na(ann_map$GENENAME)), ]
ann_map <- ann_map[!duplicated(ann_map$ENSEMBL), ]
annotation$SYMBOL <- ann_map$SYMBOL[match(annotation$ensembl_id, ann_map$ENSEMBL)]
annotation$GENENAME <- ann_map$GENENAME[match(annotation$ensembl_id, ann_map$ENSEMBL)]
annotation$ENTREZID <- ann_map$ENTREZID[match(annotation$ensembl_id, ann_map$ENSEMBL)]
annotation$gene_label <- ifelse(!is.na(annotation$SYMBOL) & annotation$SYMBOL != "",
                                annotation$SYMBOL, annotation$ensembl_id)
write.csv(annotation, file.path(table_dir, "gene_annotation.csv"), row.names = FALSE)

message("Preparing primary limma-voom model.")
sample_meta_model <- sample_meta
sample_meta_model$condition <- relevel(sample_meta_model$condition, ref = "C")
design <- model.matrix(~ pcs + condition, data = sample_meta_model)
rownames(design) <- sample_meta_model$sample_id

nonzero <- rowSums(agg_counts) > 0
y <- DGEList(
  counts = agg_counts[nonzero, sample_meta_model$sample_id, drop = FALSE],
  genes = annotation[nonzero, , drop = FALSE]
)
keep <- filterByExpr(y, design = design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")

filtering_summary <- data.frame(
  metric = c("input_genes", "nonzero_genes", "genes_after_filterByExpr", "biological_samples", "sequencing_lanes"),
  value = c(nrow(agg_counts), sum(nonzero), sum(keep), ncol(agg_counts), ncol(lane_counts))
)
write.csv(filtering_summary, file.path(table_dir, "filtering_summary.csv"), row.names = FALSE)

png(file.path(figure_dir, "qc_voom_mean_variance.png"), width = 7, height = 5, units = "in", res = 300)
voom_obj <- voomWithQualityWeights(y, design = design, plot = TRUE)
dev.off()
pdf(file.path(figure_dir, "qc_voom_mean_variance.pdf"), width = 7, height = 5)
invisible(voomWithQualityWeights(y, design = design, plot = TRUE))
dev.off()

fit <- lmFit(voom_obj, design)
fit <- eBayes(fit, robust = TRUE)

voom_logcpm <- voom_obj$E
write.csv(voom_logcpm, file.path(table_dir, "tmm_voom_logCPM_filtered_genes.csv"), row.names = TRUE)

base_theme <- theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold")
  )

condition_colors <- c(C = "#4D4D4D", `2DG` = "#D9822B", LDHA = "#2F6F9F")
deg_colors <- c(Down = "#2F6F9F", NS = "#BDBDBD", Up = "#D9822B", `FDR only` = "#7A7A7A")

save_gg <- function(plot, stem, width = 7, height = 5) {
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = width, height = height,
         units = "in", dpi = 300, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height,
         units = "in", device = cairo_pdf, bg = "white")
}

message("Generating QC figures.")
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
  scale_fill_manual(values = condition_colors, labels = c(C = "Control", `2DG` = "2DG + OX", LDHA = "LDHA + OX")) +
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
save_gg(p_lib, "qc_library_sizes", width = 8, height = 4.8)

gene_vars <- apply(voom_logcpm, 1, var)
pca_genes <- names(sort(gene_vars, decreasing = TRUE))[seq_len(min(5000, length(gene_vars)))]
pca <- prcomp(t(voom_logcpm[pca_genes, , drop = FALSE]), center = TRUE, scale. = FALSE)
pve <- (pca$sdev^2) / sum(pca$sdev^2)
pca_df <- cbind(sample_meta_model, as.data.frame(pca$x[, 1:4, drop = FALSE]))
  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, shape = pcs, label = sample_id)) +
  geom_hline(yintercept = 0, color = "grey88", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey88", linewidth = 0.3) +
  geom_point(size = 3.3, stroke = 0.7) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = Inf, box.padding = 0.35, min.segment.length = 0,
                           show.legend = FALSE) +
  scale_color_manual(values = condition_colors, labels = c(C = "Control", `2DG` = "2DG + OX", LDHA = "LDHA + OX")) +
  labs(
    title = "PCA of TMM-normalized voom expression",
    subtitle = "Top 5,000 most variable retained genes; paired PCS identity is shown by point shape",
    x = paste0("PC1 (", round(100 * pve[1], 1), "%)"),
    y = paste0("PC2 (", round(100 * pve[2], 1), "%)"),
    color = "Condition",
    shape = "PCS block"
  ) +
  base_theme
save_gg(p_pca, "qc_pca_voom_top_variable_genes", width = 7.4, height = 5.8)

sample_cor <- cor(voom_logcpm[pca_genes, , drop = FALSE], method = "pearson")
ann_col <- data.frame(
  PCS = sample_meta_model$pcs,
  Condition = sample_meta_model$condition
)
rownames(ann_col) <- sample_meta_model$sample_id
ann_colors <- list(
  Condition = condition_colors,
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
  main = "Sample correlation: voom logCPM",
  filename = file.path(figure_dir, "qc_sample_correlation_heatmap.png"),
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
  main = "Sample correlation: voom logCPM",
  filename = file.path(figure_dir, "qc_sample_correlation_heatmap.pdf"),
  width = 7.5,
  height = 6.8
)

message("Extracting differential expression results.")
contrast_specs <- list(
  `2DG_vs_C` = list(coef = "condition2DG", label = "2DG + OX vs Control", condition = "2DG"),
  `LDHA_vs_C` = list(coef = "conditionLDHA", label = "LDHA + OX vs Control", condition = "LDHA")
)

condition_means <- data.frame(
  gene_id = rownames(agg_counts),
  mean_count_C = rowMeans(agg_counts[, sample_meta_model$condition == "C", drop = FALSE]),
  mean_count_2DG = rowMeans(agg_counts[, sample_meta_model$condition == "2DG", drop = FALSE]),
  mean_count_LDHA = rowMeans(agg_counts[, sample_meta_model$condition == "LDHA", drop = FALSE]),
  stringsAsFactors = FALSE
)

classify_deg <- function(res) {
  sig <- rep("NS", nrow(res))
  sig[res$adj.P.Val < 0.05 & abs(res$logFC) < 1] <- "FDR only"
  sig[res$adj.P.Val < 0.05 & res$logFC >= 1] <- "Up"
  sig[res$adj.P.Val < 0.05 & res$logFC <= -1] <- "Down"
  factor(sig, levels = c("Down", "NS", "FDR only", "Up"))
}

limma_results <- list()
for (nm in names(contrast_specs)) {
  spec <- contrast_specs[[nm]]
  tt <- topTable(fit, coef = spec$coef, number = Inf, sort.by = "P")
  tt$gene_id <- rownames(tt)
  carried_annotation <- intersect(colnames(tt), setdiff(colnames(annotation), "gene_id"))
  if (length(carried_annotation) > 0) {
    tt <- tt[, setdiff(colnames(tt), carried_annotation), drop = FALSE]
  }
  tt <- merge(tt, annotation, by = "gene_id", all.x = TRUE, sort = FALSE)
  tt <- merge(tt, condition_means, by = "gene_id", all.x = TRUE, sort = FALSE)
  tt <- tt[match(rownames(topTable(fit, coef = spec$coef, number = Inf, sort.by = "P")), tt$gene_id), ]
  tt$contrast <- nm
  tt$contrast_label <- spec$label
  tt$direction <- classify_deg(tt)
  tt$neg_log10_fdr <- -log10(pmax(tt$adj.P.Val, .Machine$double.xmin))
  tt$rank_metric_t <- tt$t
  limma_results[[nm]] <- tt
  write.csv(tt, file.path(table_dir, paste0("limma_voom_", nm, "_all_genes.csv")), row.names = FALSE)
  write.csv(tt[tt$adj.P.Val < 0.05 & abs(tt$logFC) >= 1, ],
            file.path(table_dir, paste0("limma_voom_", nm, "_DEG_FDR05_abslog2FC1.csv")), row.names = FALSE)
}

deg_summary <- do.call(rbind, lapply(names(limma_results), function(nm) {
  res <- limma_results[[nm]]
  data.frame(
    contrast = nm,
    contrast_label = contrast_specs[[nm]]$label,
    genes_tested = nrow(res),
    fdr_0_05 = sum(res$adj.P.Val < 0.05, na.rm = TRUE),
    up_fdr_0_05_abslog2fc_1 = sum(res$adj.P.Val < 0.05 & res$logFC >= 1, na.rm = TRUE),
    down_fdr_0_05_abslog2fc_1 = sum(res$adj.P.Val < 0.05 & res$logFC <= -1, na.rm = TRUE),
    deg_fdr_0_05_abslog2fc_1 = sum(res$adj.P.Val < 0.05 & abs(res$logFC) >= 1, na.rm = TRUE)
  )
}))
write.csv(deg_summary, file.path(table_dir, "limma_voom_contrast_summary.csv"), row.names = FALSE)

plot_volcano <- function(res, stem, title) {
  labels <- res[res$adj.P.Val < 0.05 & abs(res$logFC) >= 1, ]
  labels <- labels[order(labels$adj.P.Val, -abs(labels$logFC)), ]
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
      subtitle = "Primary paired limma-voom model; dashed lines mark FDR 0.05 and |log2FC| 1",
      x = "log2 fold-change",
      y = "-log10 adjusted P value",
      color = "Class"
    ) +
    base_theme
  save_gg(p, stem, width = 7, height = 5.8)
}

plot_ma <- function(res, stem, title) {
  p <- ggplot(res, aes(x = AveExpr, y = logFC, color = direction)) +
    geom_hline(yintercept = 0, color = "grey45", linewidth = 0.35) +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "grey65", linewidth = 0.3) +
    geom_point(size = 1.0, alpha = 0.72, stroke = 0) +
    scale_color_manual(values = deg_colors, drop = TRUE) +
    labs(
      title = title,
      subtitle = "Average expression is voom average logCPM; dashed lines mark |log2FC| 1",
      x = "Average logCPM",
      y = "log2 fold-change",
      color = "Class"
    ) +
    base_theme
  save_gg(p, stem, width = 7, height = 5.2)
}

for (nm in names(limma_results)) {
  plot_volcano(limma_results[[nm]], paste0("deg_volcano_", nm), paste0("Volcano plot: ", contrast_specs[[nm]]$label))
  plot_ma(limma_results[[nm]], paste0("deg_ma_", nm), paste0("MA plot: ", contrast_specs[[nm]]$label))
}

deg_long <- do.call(rbind, lapply(names(limma_results), function(nm) {
  res <- limma_results[[nm]]
  data.frame(contrast = contrast_specs[[nm]]$label, direction = as.character(res$direction), stringsAsFactors = FALSE)
}))
deg_bar_df <- as.data.frame(table(deg_long$contrast, deg_long$direction), stringsAsFactors = FALSE)
colnames(deg_bar_df) <- c("contrast", "direction", "n")
deg_bar_df <- deg_bar_df[deg_bar_df$direction %in% c("Up", "Down", "FDR only"), ]
p_deg_bar <- ggplot(deg_bar_df, aes(x = contrast, y = n, fill = direction)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65, color = "grey25", linewidth = 0.2) +
  geom_text(aes(label = n), position = position_dodge(width = 0.75), vjust = -0.35, size = 3) +
  scale_fill_manual(values = deg_colors[c("Down", "FDR only", "Up")], drop = FALSE) +
  labs(
    title = "Differential expression burden by treatment",
    subtitle = "Up/down require FDR < 0.05 and |log2FC| >= 1; FDR only has smaller absolute effect",
    x = NULL,
    y = "Number of genes",
    fill = "Class"
  ) +
  coord_cartesian(ylim = c(0, max(deg_bar_df$n) * 1.15 + 1)) +
  base_theme
save_gg(p_deg_bar, "deg_summary_bar", width = 7, height = 4.8)

message("Generating DEG heatmaps.")
top_gene_ids <- unique(unlist(lapply(limma_results, function(res) {
  sig <- res[res$adj.P.Val < 0.05 & abs(res$logFC) >= 1, ]
  if (nrow(sig) < 20) sig <- res
  head(sig$gene_id[order(sig$adj.P.Val, -abs(sig$logFC))], 50)
})))
top_gene_ids <- intersect(top_gene_ids, rownames(voom_logcpm))
heat_mat <- voom_logcpm[top_gene_ids, sample_meta_model$sample_id, drop = FALSE]
heat_labels <- annotation$gene_label[match(rownames(heat_mat), annotation$gene_id)]
heat_labels <- make.unique(ifelse(is.na(heat_labels) | heat_labels == "", rownames(heat_mat), heat_labels))
rownames(heat_mat) <- heat_labels
ann_col_heat <- data.frame(
  PCS = sample_meta_model$pcs,
  Condition = sample_meta_model$condition
)
rownames(ann_col_heat) <- sample_meta_model$sample_id
pheatmap(
  heat_mat,
  scale = "row",
  annotation_col = ann_col_heat,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("#2166ac", "#f7f7f7", "#b35806"))(100),
  border_color = NA,
  fontsize_row = ifelse(nrow(heat_mat) > 80, 5, 7),
  fontsize_col = 8,
  main = "Top differential genes across both contrasts",
  filename = file.path(figure_dir, "deg_top_genes_heatmap.png"),
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
  main = "Top differential genes across both contrasts",
  filename = file.path(figure_dir, "deg_top_genes_heatmap.pdf"),
  width = 8.5,
  height = max(6.5, min(14, 2 + nrow(heat_mat) * 0.12))
)

message("Running DESeq2 sensitivity analysis on rounded lane-summed counts.")
deseq_results <- list()
deseq_summary <- data.frame()
deseq_error <- NULL
tryCatch({
  rounded_counts <- round(agg_counts[rownames(y), sample_meta_model$sample_id, drop = FALSE])
  dds <- DESeqDataSetFromMatrix(
    countData = rounded_counts,
    colData = sample_meta_model,
    design = ~ pcs + condition
  )
  dds <- DESeq(dds, quiet = TRUE)
  for (nm in names(contrast_specs)) {
    cond <- contrast_specs[[nm]]$condition
    dr <- as.data.frame(results(dds, contrast = c("condition", cond, "C"), alpha = 0.05))
    dr$gene_id <- rownames(dr)
    dr <- merge(dr, annotation, by = "gene_id", all.x = TRUE, sort = FALSE)
    dr <- dr[match(rownames(results(dds, contrast = c("condition", cond, "C"), alpha = 0.05)), dr$gene_id), ]
    dr$contrast <- nm
    dr$contrast_label <- contrast_specs[[nm]]$label
    dr$direction <- factor(ifelse(!is.na(dr$padj) & dr$padj < 0.05 & dr$log2FoldChange >= 1, "Up",
                                  ifelse(!is.na(dr$padj) & dr$padj < 0.05 & dr$log2FoldChange <= -1, "Down",
                                         ifelse(!is.na(dr$padj) & dr$padj < 0.05, "FDR only", "NS"))),
                           levels = c("Down", "NS", "FDR only", "Up"))
    deseq_results[[nm]] <- dr
    write.csv(dr, file.path(table_dir, paste0("deseq2_rounded_counts_", nm, "_all_genes.csv")), row.names = FALSE)
  }
  deseq_summary <- do.call(rbind, lapply(names(deseq_results), function(nm) {
    dr <- deseq_results[[nm]]
    lr <- limma_results[[nm]]
    merged <- merge(lr[, c("gene_id", "logFC")], dr[, c("gene_id", "log2FoldChange")], by = "gene_id")
    data.frame(
      contrast = nm,
      contrast_label = contrast_specs[[nm]]$label,
      deseq2_fdr_0_05 = sum(!is.na(dr$padj) & dr$padj < 0.05),
      deseq2_deg_fdr_0_05_abslog2fc_1 = sum(!is.na(dr$padj) & dr$padj < 0.05 & abs(dr$log2FoldChange) >= 1),
      logfc_spearman_vs_limma = suppressWarnings(cor(merged$logFC, merged$log2FoldChange, method = "spearman", use = "complete.obs")),
      logfc_pearson_vs_limma = suppressWarnings(cor(merged$logFC, merged$log2FoldChange, method = "pearson", use = "complete.obs"))
    )
  }))
  write.csv(deseq_summary, file.path(table_dir, "deseq2_sensitivity_summary.csv"), row.names = FALSE)

  for (nm in names(deseq_results)) {
    lr <- limma_results[[nm]][, c("gene_id", "gene_label", "logFC", "adj.P.Val")]
    dr <- deseq_results[[nm]][, c("gene_id", "log2FoldChange", "padj")]
    cmp <- merge(lr, dr, by = "gene_id")
    cmp <- cmp[is.finite(cmp$logFC) & is.finite(cmp$log2FoldChange), ]
    rho <- suppressWarnings(cor(cmp$logFC, cmp$log2FoldChange, method = "spearman", use = "complete.obs"))
    p <- ggplot(cmp, aes(x = logFC, y = log2FoldChange)) +
      geom_hline(yintercept = 0, color = "grey88", linewidth = 0.3) +
      geom_vline(xintercept = 0, color = "grey88", linewidth = 0.3) +
      geom_point(size = 0.8, alpha = 0.45, color = "#2F6F9F", stroke = 0) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.35) +
      labs(
        title = paste0("Sensitivity: DESeq2 vs limma-voom log2FC (", contrast_specs[[nm]]$label, ")"),
        subtitle = paste0("DESeq2 used rounded lane-summed counts; Spearman rho = ", round(rho, 3)),
        x = "limma-voom log2FC",
        y = "DESeq2 log2FC"
      ) +
      base_theme
    save_gg(p, paste0("sensitivity_deseq2_vs_limma_logfc_", nm), width = 5.8, height = 5.4)
  }
}, error = function(e) {
  deseq_error <<- conditionMessage(e)
  writeLines(deseq_error, file.path(log_dir, "deseq2_sensitivity_error.txt"))
})

message("Running GO enrichment and GSEA.")
make_ranked_entrez <- function(res) {
  x <- res[!is.na(res$ENTREZID) & is.finite(res$t), c("ENTREZID", "t")]
  x <- x[order(abs(x$t), decreasing = TRUE), ]
  x <- x[!duplicated(x$ENTREZID), ]
  ranks <- x$t
  names(ranks) <- x$ENTREZID
  sort(ranks, decreasing = TRUE)
}

gsea_results <- list()
ora_results <- list()
universe_entrez <- unique(na.omit(annotation$ENTREZID[match(rownames(y), annotation$gene_id)]))
for (nm in names(limma_results)) {
  res <- limma_results[[nm]]
  ranks <- make_ranked_entrez(res)
  gsea <- tryCatch({
    suppressMessages(
      gseGO(
        geneList = ranks,
        OrgDb = org.Hs.eg.db,
        ont = "BP",
        keyType = "ENTREZID",
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        verbose = FALSE,
        eps = 0
      )
    )
  }, error = function(e) {
    writeLines(conditionMessage(e), file.path(log_dir, paste0("gseGO_error_", nm, ".txt")))
    NULL
  })
  if (!is.null(gsea)) {
    gsea_df <- as.data.frame(gsea)
    gsea_results[[nm]] <- gsea_df
    write.csv(gsea_df, file.path(table_dir, paste0("go_bp_gsea_", nm, ".csv")), row.names = FALSE)
  } else {
    gsea_results[[nm]] <- data.frame()
  }

  sig_up <- unique(na.omit(res$ENTREZID[res$adj.P.Val < 0.05 & res$logFC >= 1]))
  sig_down <- unique(na.omit(res$ENTREZID[res$adj.P.Val < 0.05 & res$logFC <= -1]))
  for (direction in c("up", "down")) {
    genes <- if (direction == "up") sig_up else sig_down
    ora_nm <- paste(nm, direction, sep = "_")
    if (length(genes) >= 10) {
      ora <- tryCatch({
        suppressMessages(
          enrichGO(
            gene = genes,
            universe = universe_entrez,
            OrgDb = org.Hs.eg.db,
            keyType = "ENTREZID",
            ont = "BP",
            pAdjustMethod = "BH",
            pvalueCutoff = 1,
            qvalueCutoff = 1,
            minGSSize = 10,
            maxGSSize = 500,
            readable = TRUE
          )
        )
      }, error = function(e) {
        writeLines(conditionMessage(e), file.path(log_dir, paste0("enrichGO_error_", ora_nm, ".txt")))
        NULL
      })
      ora_df <- if (!is.null(ora)) as.data.frame(ora) else data.frame()
    } else {
      ora_df <- data.frame()
    }
    ora_results[[ora_nm]] <- ora_df
    write.csv(ora_df, file.path(table_dir, paste0("go_bp_ora_", ora_nm, ".csv")), row.names = FALSE)
  }
}

plot_gsea <- function(df, stem, title) {
  if (nrow(df) == 0) {
    p <- ggplot() +
      annotate("text", x = 0, y = 0, label = "No GO BP GSEA results available", size = 4) +
      xlim(-1, 1) + ylim(-1, 1) +
      labs(title = title, x = NULL, y = NULL) +
      theme_void(base_size = 11)
    save_gg(p, stem, width = 7, height = 4.5)
    return(invisible(NULL))
  }
  df <- df[is.finite(df$NES) & !is.na(df$p.adjust), ]
  pos <- head(df[df$NES > 0, ][order(df$p.adjust[df$NES > 0], -abs(df$NES[df$NES > 0])), ], 8)
  neg <- head(df[df$NES < 0, ][order(df$p.adjust[df$NES < 0], -abs(df$NES[df$NES < 0])), ], 8)
  top <- rbind(pos, neg)
  if (nrow(top) == 0) top <- head(df[order(df$p.adjust), ], 16)
  top$Description <- factor(top$Description, levels = top$Description[order(top$NES)])
  top$direction <- ifelse(top$NES > 0, "Enriched in treatment", "Enriched in control")
  p <- ggplot(top, aes(x = NES, y = Description, size = setSize, color = p.adjust)) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.35) +
    geom_point(alpha = 0.9) +
    scale_color_gradient(low = "#D9822B", high = "#2F6F9F", trans = "reverse") +
    labs(
      title = title,
      subtitle = "Top positive and negative GO Biological Process GSEA terms",
      x = "Normalized enrichment score",
      y = NULL,
      color = "Adjusted P",
      size = "Gene set size"
    ) +
    base_theme
  save_gg(p, stem, width = 8.5, height = 6.3)
}

for (nm in names(gsea_results)) {
  plot_gsea(gsea_results[[nm]], paste0("pathway_go_bp_gsea_", nm),
            paste0("GO BP GSEA: ", contrast_specs[[nm]]$label))
}

message("Writing Excel workbook.")
wb <- createWorkbook()
addWorksheet(wb, "sample_metadata")
writeDataTable(wb, "sample_metadata", sample_meta)
addWorksheet(wb, "lane_metadata")
writeDataTable(wb, "lane_metadata", lane_meta)
addWorksheet(wb, "filtering_summary")
writeDataTable(wb, "filtering_summary", filtering_summary)
addWorksheet(wb, "limma_summary")
writeDataTable(wb, "limma_summary", deg_summary)
if (nrow(deseq_summary) > 0) {
  addWorksheet(wb, "deseq2_sensitivity")
  writeDataTable(wb, "deseq2_sensitivity", deseq_summary)
}
for (nm in names(limma_results)) {
  addWorksheet(wb, paste0("limma_", substr(nm, 1, 24)))
  writeDataTable(wb, paste0("limma_", substr(nm, 1, 24)), limma_results[[nm]])
  sig_sheet <- paste0("DEG_", substr(nm, 1, 25))
  addWorksheet(wb, sig_sheet)
  writeDataTable(wb, sig_sheet, limma_results[[nm]][limma_results[[nm]]$adj.P.Val < 0.05 & abs(limma_results[[nm]]$logFC) >= 1, ])
  gsea_sheet <- paste0("GSEA_", substr(nm, 1, 24))
  addWorksheet(wb, gsea_sheet)
  writeDataTable(wb, gsea_sheet, gsea_results[[nm]])
}
saveWorkbook(wb, file.path(table_dir, "bulk_rnaseq_2dg_ldha_publication_results.xlsx"), overwrite = TRUE)

top_genes_summary <- do.call(rbind, lapply(names(limma_results), function(nm) {
  res <- limma_results[[nm]]
  sig <- res[res$adj.P.Val < 0.05 & abs(res$logFC) >= 1, ]
  if (nrow(sig) == 0) sig <- head(res[order(res$adj.P.Val), ], 10)
  sig <- sig[order(sig$adj.P.Val, -abs(sig$logFC)), ]
  head(sig[, c("contrast_label", "gene_id", "gene_label", "SYMBOL", "GENENAME", "logFC", "AveExpr", "P.Value", "adj.P.Val", "direction")], 12)
}))
write.csv(top_genes_summary, file.path(table_dir, "top_genes_summary_for_report.csv"), row.names = FALSE)

top_gsea_summary <- do.call(rbind, lapply(names(gsea_results), function(nm) {
  df <- gsea_results[[nm]]
  if (nrow(df) == 0) return(data.frame())
  df <- df[order(df$p.adjust, -abs(df$NES)), ]
  df <- head(df[, c("ID", "Description", "setSize", "NES", "pvalue", "p.adjust", "qvalue")], 10)
  df$contrast <- nm
  df$contrast_label <- contrast_specs[[nm]]$label
  df
}))
write.csv(top_gsea_summary, file.path(table_dir, "top_go_gsea_summary_for_report.csv"), row.names = FALSE)

source_notes <- list(
  timestamp = timestamp,
  input_file = input_file,
  analysis_dir = analysis_dir,
  primary_model = "edgeR TMM normalization + limma voomWithQualityWeights + paired design (~ PCS + condition), robust empirical Bayes",
  sensitivity_model = if (is.null(deseq_error)) "DESeq2 paired design on rounded technical-lane-summed counts" else paste("DESeq2 sensitivity failed:", deseq_error),
  differential_expression_threshold = "FDR < 0.05 and |log2FC| >= 1 for high-confidence DEG calls; FDR-only genes are also reported",
  pathway_model = "clusterProfiler GO Biological Process GSEA on moderated t-statistics and enrichGO ORA for up/down DEG sets",
  caveats = c(
    "The provided matrix contains fractional counts, consistent with estimated gene-level counts. The primary model therefore uses limma-voom, which accepts non-integer abundance estimates.",
    "Technical lanes were summed before modeling to avoid treating lanes as biological replicates.",
    "No FASTQ/BAM alignment QC metrics, gene length offsets, or batch annotations beyond PCS block were provided."
  )
)
saveRDS(
  list(
    source_notes = source_notes,
    lane_meta = lane_meta,
    sample_meta = sample_meta,
    filtering_summary = filtering_summary,
    deg_summary = deg_summary,
    deseq_summary = deseq_summary,
    top_genes_summary = top_genes_summary,
    top_gsea_summary = top_gsea_summary,
    limma_results = limma_results,
    gsea_results = gsea_results,
    pca_variance = pve,
    design = design
  ),
  file.path(object_dir, "analysis_summary.rds")
)

writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo.txt"))

message("Writing and rendering technical report.")
rmd_path <- file.path(report_dir, "bulk_rnaseq_2dg_ldha_report.Rmd")
report_lines <- c(
  "---",
  "title: \"Bulk RNA-seq differential expression: 2DG + OX and LDHA + OX versus control\"",
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
  "summary <- readRDS(file.path(out_dir, 'objects', 'analysis_summary.rds'))",
  "fmt_int <- function(x) format(x, big.mark = ',', scientific = FALSE, trim = TRUE)",
  "fmt_num <- function(x, digits = 3) format(round(x, digits), nsmall = digits, trim = TRUE)",
  "fig <- function(name) knitr::include_graphics(file.path(out_dir, 'figures', name))",
  "deg_summary <- summary$deg_summary",
  "top_genes <- summary$top_genes_summary",
  "top_gsea <- summary$top_gsea_summary",
  "filtering <- summary$filtering_summary",
  "sample_meta <- summary$sample_meta",
  "deseq_summary <- summary$deseq_summary",
  "```",
  "",
  "## Technical Summary",
  "",
  "This analysis compares two treatment conditions against paired controls: `2DG` is interpreted as `2DG + OX`, and `LDHA` is interpreted as `LDHA + OX`. The input contained `r fmt_int(filtering$value[filtering$metric == 'input_genes'])` gene rows across `r fmt_int(filtering$value[filtering$metric == 'sequencing_lanes'])` sequencing-lane columns. Technical lanes belonging to the same PCS-condition sample were summed, yielding `r fmt_int(nrow(sample_meta))` biological profiles: five PCS blocks measured under control, 2DG + OX, and LDHA + OX.",
  "",
  "The primary differential-expression model uses TMM normalization and `limma-voom` with sample quality weights in a paired/block design: `~ PCS + condition`. This is the most defensible primary model for this matrix because counts are fractional rather than strict raw integers. High-confidence DEG calls are defined as FDR < 0.05 and absolute log2 fold-change >= 1; genes passing FDR alone are retained in the all-gene tables.",
  "",
  "```{r summary-table}",
  "show_summary <- deg_summary[, c('contrast_label', 'genes_tested', 'fdr_0_05', 'up_fdr_0_05_abslog2fc_1', 'down_fdr_0_05_abslog2fc_1', 'deg_fdr_0_05_abslog2fc_1')]",
  "colnames(show_summary) <- c('Contrast', 'Genes tested', 'FDR < 0.05', 'Up DEG', 'Down DEG', 'Total DEG')",
  "knitr::kable(show_summary, digits = 3)",
  "```",
  "",
  "## The paired samples support treatment-level comparisons after lane consolidation",
  "",
  "The sample set is balanced at the biological level after technical lane consolidation: each PCS block has one control, one 2DG + OX, and one LDHA + OX profile. The library-size plot is therefore a QC check for depth and technical-lane consolidation rather than evidence of biological replication.",
  "",
  "```{r qc-library, out.width='100%'}",
  "fig('qc_library_sizes.png')",
  "```",
  "",
  "The PCA uses the top 5,000 most variable retained genes after TMM/voom transformation. The plot should be read as a sample-level QC view: treatment separation is interpretable only after accounting for the paired PCS structure in the model.",
  "",
  "```{r qc-pca, out.width='100%'}",
  "fig('qc_pca_voom_top_variable_genes.png')",
  "```",
  "",
  "The correlation heatmap checks whether any sample behaves as a strong outlier after normalization. Correlation structure is expected to reflect both paired PCS identity and treatment response.",
  "",
  "```{r qc-cor, out.width='100%'}",
  "fig('qc_sample_correlation_heatmap.png')",
  "```",
  "",
  "## Both treatment contrasts were tested with the same PCS-blocked model",
  "",
  "The DEG summary below separates high-confidence up/down genes from genes that meet FDR but have smaller absolute effect sizes. This avoids mixing statistical detectability with a publication-scale effect threshold.",
  "",
  "```{r deg-summary-plot, out.width='100%'}",
  "fig('deg_summary_bar.png')",
  "```",
  "",
  "### 2DG + OX versus control",
  "",
  "The volcano and MA plots show the same contrast from complementary angles: the volcano emphasizes significance and effect size, while the MA plot checks whether fold-change estimates depend strongly on expression level.",
  "",
  "```{r volcano-2dg, out.width='100%'}",
  "fig('deg_volcano_2DG_vs_C.png')",
  "```",
  "",
  "```{r ma-2dg, out.width='100%'}",
  "fig('deg_ma_2DG_vs_C.png')",
  "```",
  "",
  "### LDHA + OX versus control",
  "",
  "The LDHA + OX contrast was modeled with the identical normalization, filter, and paired design, so effect sizes and FDR values are directly comparable to the 2DG + OX contrast.",
  "",
  "```{r volcano-ldha, out.width='100%'}",
  "fig('deg_volcano_LDHA_vs_C.png')",
  "```",
  "",
  "```{r ma-ldha, out.width='100%'}",
  "fig('deg_ma_LDHA_vs_C.png')",
  "```",
  "",
  "## Top differential genes provide a manuscript-ready shortlist",
  "",
  "The heatmap shows the strongest differential genes across both contrasts, using row-scaled voom logCPM values. It is intended as a compact figure-panel candidate; exact statistics should be taken from the exported all-gene tables.",
  "",
  "```{r heatmap, out.width='100%'}",
  "fig('deg_top_genes_heatmap.png')",
  "```",
  "",
  "```{r top-gene-table}",
  "tg <- top_genes[, c('contrast_label', 'gene_label', 'GENENAME', 'logFC', 'AveExpr', 'adj.P.Val', 'direction')]",
  "colnames(tg) <- c('Contrast', 'Gene', 'Gene name', 'log2FC', 'Ave logCPM', 'FDR', 'Class')",
  "knitr::kable(tg, digits = 4)",
  "```",
  "",
  "## GO Biological Process analysis summarizes pathway-level signal",
  "",
  "GO Biological Process GSEA was run on the full ranked gene list for each contrast, using moderated t-statistics as the ranking metric. Positive NES values indicate enrichment toward the treatment-up end of the ranked list; negative NES values indicate enrichment toward control.",
  "",
  "```{r gsea-2dg, out.width='100%'}",
  "fig('pathway_go_bp_gsea_2DG_vs_C.png')",
  "```",
  "",
  "```{r gsea-ldha, out.width='100%'}",
  "fig('pathway_go_bp_gsea_LDHA_vs_C.png')",
  "```",
  "",
  "```{r gsea-table}",
  "if (nrow(top_gsea) > 0) {",
  "  gsea_show <- top_gsea[, c('contrast_label', 'Description', 'setSize', 'NES', 'p.adjust', 'qvalue')]",
  "  colnames(gsea_show) <- c('Contrast', 'GO Biological Process', 'Set size', 'NES', 'FDR', 'q value')",
  "  knitr::kable(gsea_show, digits = 4)",
  "} else {",
  "  cat('No GO GSEA results were returned.')",
  "}",
  "```",
  "",
  "## Scope, data, and definitions",
  "",
  "- Biological comparison: 2DG + OX vs control and LDHA + OX vs control.",
  "- Biological blocking factor: PCS sample identity, modeled as a fixed paired effect.",
  "- Technical lanes: sequencing-lane columns with the same PCS and condition were summed before modeling.",
  "- Primary expression unit: TMM-normalized voom logCPM for QC, heatmaps, and linear modeling.",
  "- DEG threshold: FDR < 0.05 and |log2FC| >= 1 for the publication shortlist; all tested genes are exported without effect-size filtering.",
  "- Gene identifiers: Ensembl IDs were version-stripped for annotation through `org.Hs.eg.db`; unmapped genes remain in the analysis with their Ensembl IDs.",
  "",
  "```{r sample-table}",
  "sm <- sample_meta[, c('sample_id', 'pcs', 'condition_label', 'n_technical_lanes', 'library_size')]",
  "colnames(sm) <- c('Sample', 'PCS block', 'Condition', 'Technical lanes', 'Summed library size')",
  "knitr::kable(sm)",
  "```",
  "",
  "## Methodology",
  "",
  "1. Imported the Galaxy gene-level matrix and parsed sample names into PCS block, condition, and technical-lane metadata.",
  "2. Summed technical lanes belonging to the same PCS-condition biological sample.",
  "3. Removed genes with zero total signal, then applied `edgeR::filterByExpr` using the paired design matrix.",
  "4. Normalized retained genes with TMM and fitted `limma::voomWithQualityWeights` using `~ PCS + condition`.",
  "5. Tested `condition2DG` and `conditionLDHA` coefficients against control with robust empirical Bayes moderation.",
  "6. Annotated genes with `org.Hs.eg.db`, exported all-gene and high-confidence DEG tables, and generated publication-style QC and DEG figures.",
  "7. Ran a DESeq2 sensitivity analysis on rounded lane-summed counts. This is secondary because the source matrix contains fractional counts.",
  "8. Ran GO Biological Process GSEA with `clusterProfiler::gseGO` on moderated t-statistics, plus GO over-representation tables for up/down DEG sets when enough genes were available.",
  "",
  "```{r deseq-sensitivity}",
  "if (nrow(deseq_summary) > 0) {",
  "  ds <- deseq_summary[, c('contrast_label', 'deseq2_fdr_0_05', 'deseq2_deg_fdr_0_05_abslog2fc_1', 'logfc_spearman_vs_limma', 'logfc_pearson_vs_limma')]",
  "  colnames(ds) <- c('Contrast', 'DESeq2 FDR < 0.05', 'DESeq2 DEG', 'Spearman log2FC vs limma', 'Pearson log2FC vs limma')",
  "  knitr::kable(ds, digits = 4)",
  "} else {",
  "  cat(summary$source_notes$sensitivity_model)",
  "}",
  "```",
  "",
  "## Limitations, uncertainty, and robustness checks",
  "",
  "- The supplied matrix contains fractional count values. The primary limma-voom model is appropriate for non-integer estimated counts, while DESeq2 is provided only as a rounded-count sensitivity analysis.",
  "- No raw sequencing QC, alignment metrics, strandedness checks, or gene-length offsets were provided. The report therefore validates the count matrix and modeled contrasts, not upstream read processing.",
  "- PCS is modeled as a fixed paired effect. With five biological blocks per condition, large effects are more reliable than borderline FDR-only calls.",
  "- GO terms can be redundant and broad. The exported GSEA and ORA tables should be used to choose a concise pathway narrative for the manuscript.",
  "",
  "## Recommended next steps",
  "",
  "1. Review the DEG shortlists against expected biology and remove genes that are unannotated or low-confidence for the planned figure panel.",
  "2. If FASTQ/BAM QC is available, add upstream metrics such as mapping rate, duplication, strandedness, and gene-body coverage before submission.",
  "3. For manuscript pathway claims, collapse redundant GO terms into a small set of biological themes and validate those against the leading-edge genes.",
  "4. Confirm whether `2DG + OX` and `LDHA + OX` should be interpreted as paired perturbations in the same PCS system for all downstream text.",
  "",
  "## Further questions",
  "",
  "- Are there additional covariates such as sequencing batch, library prep batch, or RNA integrity that should be included in the design?",
  "- Should the final manuscript emphasize high-confidence |log2FC| >= 1 genes, or all FDR-significant genes regardless of effect size?",
  "- Are there curated gene sets specific to the biological system that should supplement the GO analysis?"
)
writeLines(report_lines, rmd_path)

Sys.setenv(RSTUDIO_PANDOC = project_path("external/tools"))
rmarkdown::render(
  input = rmd_path,
  output_file = "bulk_rnaseq_2dg_ldha_report.html",
  output_dir = report_dir,
  params = list(out_dir = analysis_dir),
  quiet = TRUE,
  envir = new.env(parent = globalenv())
)

message("Analysis complete.")
message("Report: ", file.path(report_dir, "bulk_rnaseq_2dg_ldha_report.html"))
message("Workbook: ", file.path(table_dir, "bulk_rnaseq_2dg_ldha_publication_results.xlsx"))
