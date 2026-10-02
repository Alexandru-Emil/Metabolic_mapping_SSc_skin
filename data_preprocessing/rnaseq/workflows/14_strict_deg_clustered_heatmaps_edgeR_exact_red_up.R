# Scientific workflow derived from Revision/RNAseq PCS met inhibitors/results/scripts/14_strict_deg_clustered_heatmaps_edgeR_exact_red_up.R
# Run from the repository root with METABOLIC_INPUT_DIR set to the external project data.
source("data_preprocessing/common/config.R")
suppressPackageStartupMessages({
  library(openxlsx)
  library(pheatmap)
})

report_file <- paste0(
  paste0(project_path("Revision/"), "/"),
  "RNAseq PCS met inhibitors/results/edgeR_batch_corrected_analysis/report/",
  "edgeR_exact_ComBatSeq_FDR05_FC15_threshold_GO_GSEA_v5/report_data.rds"
)
output_dir <- Sys.getenv(
  "RNASEQ_FIGURE_OUTPUT_DIR",
  file.path(
    getwd(),
    "output",
    "edgeR_exact_strict_DEG_heatmaps_v2_red_up_blue_down"
  )
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

report <- readRDS(report_file)
analysis_results <- report$analysis_results
gene_results <- report$gene_results
expression_long <- report$expression
sample_labels <- report$sample_labels

fdr_cutoff <- 0.05
absolute_fc_cutoff <- 1.5
absolute_log2fc_cutoff <- log2(absolute_fc_cutoff)
z_limit <- 2.5

condition_colors <- c(
  "Control" = "#59636E",
  "2DG + OX" = "#C83E3A",
  "LDHA + OX" = "#8F4E93"
)
pcs_colors <- c(
  "PCS1" = "#4477AA",
  "PCS2" = "#EE6677",
  "PCS3" = "#228833",
  "PCS4" = "#CCBB44",
  "PCS5" = "#AA3377"
)
direction_colors <- c(
  "Downregulated" = "#2166AC",
  "Upregulated" = "#B2182B"
)

extract_strict_degs <- function(contrast_key) {
  selected <- gene_results[
    gene_results$Contrast_key == contrast_key &
      is.finite(gene_results$FDR) &
      is.finite(gene_results$log2FC) &
      gene_results$FDR < fdr_cutoff &
      abs(gene_results$log2FC) > absolute_log2fc_cutoff,
    ,
    drop = FALSE
  ]
  selected <- selected[
    order(selected$FDR, -abs(selected$log2FC)),
    ,
    drop = FALSE
  ]

  if (!nrow(selected)) {
    stop("No strict DEGs found for ", contrast_key, ".")
  }
  if (anyNA(selected$Symbol) || any(!nzchar(selected$Symbol))) {
    stop("Blank symbols found in strict DEGs for ", contrast_key, ".")
  }
  if (anyDuplicated(toupper(selected$Symbol))) {
    stop("Duplicate symbols found in strict DEGs for ", contrast_key, ".")
  }
  if (anyDuplicated(selected$expression_key)) {
    stop("Duplicate expression keys found in strict DEGs for ", contrast_key, ".")
  }

  selected$Direction <- factor(
    ifelse(selected$log2FC < 0, "Downregulated", "Upregulated"),
    levels = c("Downregulated", "Upregulated")
  )
  selected
}

strict_feature_summary <- function(contrast_key) {
  selected <- analysis_results[
    analysis_results$Contrast_key == contrast_key &
      is.finite(analysis_results$FDR) &
      is.finite(analysis_results$log2FC) &
      analysis_results$FDR < fdr_cutoff &
      abs(analysis_results$log2FC) > absolute_log2fc_cutoff,
    ,
    drop = FALSE
  ]
  data.frame(
    Contrast_key = contrast_key,
    Called_feature_DEGs = nrow(selected),
    Feature_DEGs_with_symbol = sum(selected$has_symbol),
    Feature_DEGs_without_symbol = sum(!selected$has_symbol),
    stringsAsFactors = FALSE
  )
}

sample_keys_for_heatmap <- function(treatment_code) {
  keys <- names(sample_labels)
  control <- keys[grepl("_C$", keys)]
  treatment <- keys[grepl(paste0("_", treatment_code, "$"), keys)]
  expected_pcs <- paste0("PCS", 1:5)

  control <- control[match(expected_pcs, sub("_.*$", "", control))]
  treatment <- treatment[match(expected_pcs, sub("_.*$", "", treatment))]
  sample_keys <- c(control, treatment)

  if (length(sample_keys) != 10 || anyNA(sample_keys)) {
    stop("Expected five control and five ", treatment_code, " samples.")
  }
  sample_keys
}

make_expression_matrices <- function(stats, treatment_code) {
  sample_keys <- sample_keys_for_heatmap(treatment_code)
  selected_expression <- expression_long[
    expression_long$expression_key %in% stats$expression_key &
      expression_long$Sample_key %in% sample_keys,
    ,
    drop = FALSE
  ]

  matrix_logcpm <- matrix(
    NA_real_,
    nrow = nrow(stats),
    ncol = length(sample_keys),
    dimnames = list(stats$Symbol, sample_keys)
  )
  row_index <- match(selected_expression$expression_key, stats$expression_key)
  col_index <- match(selected_expression$Sample_key, sample_keys)
  matrix_logcpm[cbind(row_index, col_index)] <- selected_expression$logCPM

  if (anyNA(matrix_logcpm)) {
    stop("Missing expression values in the ", treatment_code, " heatmap.")
  }

  matrix_z <- t(apply(matrix_logcpm, 1, function(values) {
    value_sd <- stats::sd(values)
    if (!is.finite(value_sd) || value_sd == 0) {
      rep(0, length(values))
    } else {
      (values - mean(values)) / value_sd
    }
  }))
  rownames(matrix_z) <- rownames(matrix_logcpm)
  colnames(matrix_z) <- colnames(matrix_logcpm)

  correlation <- stats::cor(t(matrix_z), method = "pearson")
  if (anyNA(correlation)) {
    stop("Undefined row correlations in the ", treatment_code, " heatmap.")
  }
  distance_matrix <- 1 - correlation
  distance_matrix[distance_matrix < 0] <- 0
  diag(distance_matrix) <- 0
  row_tree <- stats::hclust(stats::as.dist(distance_matrix), method = "average")

  list(
    logCPM = matrix_logcpm,
    z_score = matrix_z,
    row_tree = row_tree
  )
}

make_column_annotation <- function(
    sample_keys,
    treatment_code,
    treatment_label) {
  condition <- ifelse(grepl("_C$", sample_keys), "Control", treatment_label)
  data.frame(
    Condition = factor(condition, levels = c("Control", treatment_label)),
    PCS = factor(sub("_.*$", "", sample_keys), levels = paste0("PCS", 1:5)),
    row.names = sample_keys,
    check.names = FALSE
  )
}

make_column_labels <- function(sample_keys, treatment_code) {
  short_condition <- ifelse(
    grepl("_C$", sample_keys),
    "C",
    ifelse(treatment_code == "2DG", "2DG", "LDHA")
  )
  paste(sub("_.*$", "", sample_keys), short_condition)
}

make_row_annotation <- function(stats) {
  data.frame(
    Direction = factor(
      stats$Direction,
      levels = c("Downregulated", "Upregulated")
    ),
    row.names = stats$Symbol,
    check.names = FALSE
  )
}

clip_z_matrix <- function(matrix_z) {
  clipped <- matrix_z
  clipped[clipped < -z_limit] <- -z_limit
  clipped[clipped > z_limit] <- z_limit
  clipped
}

build_heatmap <- function(
    matrices,
    stats,
    treatment_code,
    treatment_label,
    title,
    show_symbols,
    fontsize_row) {
  sample_keys <- colnames(matrices$z_score)
  annotation_col <- make_column_annotation(
    sample_keys,
    treatment_code,
    treatment_label
  )
  annotation_row <- make_row_annotation(stats)

  pheatmap(
    clip_z_matrix(matrices$z_score),
    color = grDevices::colorRampPalette(
      c("#2166AC", "#F7F7F7", "#B2182B")
    )(101),
    breaks = seq(-z_limit, z_limit, length.out = 102),
    cluster_rows = matrices$row_tree,
    cluster_cols = FALSE,
    gaps_col = 5,
    annotation_row = annotation_row,
    annotation_col = annotation_col,
    annotation_colors = list(
      Condition = condition_colors[c("Control", treatment_label)],
      PCS = pcs_colors,
      Direction = direction_colors
    ),
    annotation_names_row = FALSE,
    annotation_names_col = TRUE,
    show_rownames = show_symbols,
    show_colnames = TRUE,
    labels_col = make_column_labels(sample_keys, treatment_code),
    angle_col = 45,
    border_color = NA,
    fontsize = 8.2,
    fontsize_row = fontsize_row,
    fontsize_col = 7.2,
    main = paste0(
      title,
      "\nFDR < 0.05 and |FC| > 1.5",
      "\nPCS-corrected logCPM (row z-score)"
    ),
    legend_breaks = c(-2, -1, 0, 1, 2),
    legend_labels = c("-2", "-1", "0", "1", "2"),
    treeheight_row = if (show_symbols) 42 else 30,
    treeheight_col = 0,
    silent = TRUE
  )
}

save_pheatmap <- function(heatmap, stem, width_mm, height_mm) {
  pdf_file <- file.path(output_dir, paste0(stem, ".pdf"))
  png_file <- file.path(output_dir, paste0(stem, ".png"))

  grDevices::cairo_pdf(
    filename = pdf_file,
    width = width_mm / 25.4,
    height = height_mm / 25.4,
    bg = "white"
  )
  grid::grid.newpage()
  grid::grid.draw(heatmap$gtable)
  grDevices::dev.off()

  grDevices::png(
    filename = png_file,
    width = width_mm,
    height = height_mm,
    units = "mm",
    res = 600,
    type = "cairo",
    bg = "white"
  )
  grid::grid.newpage()
  grid::grid.draw(heatmap$gtable)
  grDevices::dev.off()

  c(pdf_file, png_file)
}

build_treatment_outputs <- function(
    contrast_key,
    treatment_code,
    treatment_label,
    file_prefix) {
  stats <- extract_strict_degs(contrast_key)
  matrices <- make_expression_matrices(stats, treatment_code)
  title <- paste0(
    treatment_label,
    " vs control: ",
    nrow(stats),
    " strict DEGs"
  )

  compact <- build_heatmap(
    matrices = matrices,
    stats = stats,
    treatment_code = treatment_code,
    treatment_label = treatment_label,
    title = title,
    show_symbols = FALSE,
    fontsize_row = 4
  )
  labeled <- build_heatmap(
    matrices = matrices,
    stats = stats,
    treatment_code = treatment_code,
    treatment_label = treatment_label,
    title = title,
    show_symbols = TRUE,
    fontsize_row = 5.5
  )

  labeled_height <- max(250, 62 + nrow(stats) * 1.72)
  files <- c(
    save_pheatmap(
      compact,
      paste0(file_prefix, "_strict_DEG_heatmap_compact"),
      width_mm = 150,
      height_mm = 155
    ),
    save_pheatmap(
      labeled,
      paste0(file_prefix, "_strict_DEG_heatmap_labeled"),
      width_mm = 185,
      height_mm = labeled_height
    )
  )

  cluster_order <- data.frame(
    Cluster_rank = seq_along(matrices$row_tree$order),
    Symbol = rownames(matrices$z_score)[matrices$row_tree$order],
    stringsAsFactors = FALSE
  )
  stats$Cluster_rank <- match(stats$Symbol, cluster_order$Symbol)

  list(
    stats = stats,
    matrices = matrices,
    cluster_order = cluster_order,
    files = files
  )
}

result_2dg <- build_treatment_outputs(
  contrast_key = "2DG_vs_C",
  treatment_code = "2DG",
  treatment_label = "2DG + OX",
  file_prefix = "2DG_OX_edgeR_exact"
)
result_ldha <- build_treatment_outputs(
  contrast_key = "LDHA_vs_C",
  treatment_code = "LDHA",
  treatment_label = "LDHA + OX",
  file_prefix = "LDHA_OX_edgeR_exact"
)

feature_summary_2dg <- strict_feature_summary("2DG_vs_C")
feature_summary_ldha <- strict_feature_summary("LDHA_vs_C")

if (
  feature_summary_2dg$Feature_DEGs_with_symbol != nrow(result_2dg$stats) ||
    feature_summary_ldha$Feature_DEGs_with_symbol != nrow(result_ldha$stats)
) {
  stop("Mapped feature counts do not match the unique symbol heatmap rows.")
}

stats_export_columns <- c(
  "Cluster_rank", "Symbol", "Direction", "Contrast", "log2FC", "Fold_change",
  "Statistic", "P_value", "FDR", "Average_expression", "expression_key"
)
stats_2dg_export <- result_2dg$stats[
  order(result_2dg$stats$Cluster_rank),
  stats_export_columns
]
stats_ldha_export <- result_ldha$stats[
  order(result_ldha$stats$Cluster_rank),
  stats_export_columns
]

matrix_to_export <- function(matrix_values, cluster_order) {
  ordered <- matrix_values[cluster_order$Symbol, , drop = FALSE]
  data.frame(
    Cluster_rank = cluster_order$Cluster_rank,
    Symbol = rownames(ordered),
    ordered,
    check.names = FALSE,
    row.names = NULL
  )
}

csv_2dg <- file.path(output_dir, "2DG_OX_edgeR_exact_strict_DEGs.csv")
csv_ldha <- file.path(output_dir, "LDHA_OX_edgeR_exact_strict_DEGs.csv")
xlsx_file <- file.path(output_dir, "edgeR_exact_strict_DEG_heatmap_values.xlsx")

write.csv(stats_2dg_export, csv_2dg, row.names = FALSE, na = "")
write.csv(stats_ldha_export, csv_ldha, row.names = FALSE, na = "")
notes_export <- data.frame(
  Item = c(
    "DEG method",
    "DEG threshold",
    "Heatmap expression values",
    "Heatmap cell colors",
    "DEG direction annotation",
    "2DG + OX feature-level DEGs called",
    "2DG + OX feature DEGs without mapped symbol",
    "2DG + OX unique symbols shown",
    "LDHA + OX feature-level DEGs called",
    "LDHA + OX feature DEGs without mapped symbol",
    "LDHA + OX unique symbols shown"
  ),
  Value = c(
    "edgeR exact test on ComBat-seq corrected counts",
    "FDR < 0.05 and absolute fold change > 1.5",
    "Row z-scores of ComBat-seq PCS-corrected logCPM",
    "Blue = lower relative expression; red = higher relative expression",
    "Blue = downregulated; red = upregulated",
    feature_summary_2dg$Called_feature_DEGs,
    feature_summary_2dg$Feature_DEGs_without_symbol,
    nrow(result_2dg$stats),
    feature_summary_ldha$Called_feature_DEGs,
    feature_summary_ldha$Feature_DEGs_without_symbol,
    nrow(result_ldha$stats)
  ),
  stringsAsFactors = FALSE
)
write.xlsx(
  list(
    `README` = notes_export,
    `2DG_DEGs` = stats_2dg_export,
    `LDHA_DEGs` = stats_ldha_export,
    `2DG_logCPM` = matrix_to_export(
      result_2dg$matrices$logCPM,
      result_2dg$cluster_order
    ),
    `2DG_row_zscore` = matrix_to_export(
      result_2dg$matrices$z_score,
      result_2dg$cluster_order
    ),
    `LDHA_logCPM` = matrix_to_export(
      result_ldha$matrices$logCPM,
      result_ldha$cluster_order
    ),
    `LDHA_row_zscore` = matrix_to_export(
      result_ldha$matrices$z_score,
      result_ldha$cluster_order
    )
  ),
  xlsx_file,
  asTable = TRUE,
  overwrite = TRUE
)

manifest_file <- file.path(output_dir, "figure_manifest.txt")
writeLines(
  c(
    "Figures: clustered strict-DEG heatmaps for 2DG + OX and LDHA + OX",
    "Differential-expression method: edgeR exact test",
    "DEG definition: FDR < 0.05 and absolute fold change > 1.5",
    "Expression display: row z-scores of ComBat-seq PCS-corrected logCPM",
    "Heatmap cell colors: blue = lower relative expression; red = higher relative expression",
    "DEG direction annotation: blue = downregulated; red = upregulated",
    "Columns: ten individual samples; five controls followed by five treatment samples",
    "Rows: genes clustered by 1 - Pearson correlation using average linkage",
    "Columns are not clustered",
    paste0(
      "2DG + OX strict DEGs: ",
      nrow(result_2dg$stats),
      " (",
      sum(result_2dg$stats$Direction == "Downregulated"),
      " down, ",
      sum(result_2dg$stats$Direction == "Upregulated"),
      " up)"
    ),
    paste0(
      "2DG + OX app reconciliation: ",
      feature_summary_2dg$Called_feature_DEGs,
      " feature-level DEGs called - ",
      feature_summary_2dg$Feature_DEGs_without_symbol,
      " without a mapped symbol = ",
      nrow(result_2dg$stats),
      " unique symbols shown"
    ),
    paste0(
      "LDHA + OX strict DEGs: ",
      nrow(result_ldha$stats),
      " (",
      sum(result_ldha$stats$Direction == "Downregulated"),
      " down, ",
      sum(result_ldha$stats$Direction == "Upregulated"),
      " up)"
    ),
    paste0(
      "LDHA + OX app reconciliation: ",
      feature_summary_ldha$Called_feature_DEGs,
      " feature-level DEGs called - ",
      feature_summary_ldha$Feature_DEGs_without_symbol,
      " without a mapped symbol = ",
      nrow(result_ldha$stats),
      " unique symbols shown"
    ),
    "Compact versions omit row symbols; labeled versions include every symbol",
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
  ),
  manifest_file
)

all_outputs <- c(
  result_2dg$files,
  result_ldha$files,
  csv_2dg,
  csv_ldha,
  xlsx_file,
  manifest_file
)
cat("Created:\n", paste(all_outputs, collapse = "\n"), "\n")
