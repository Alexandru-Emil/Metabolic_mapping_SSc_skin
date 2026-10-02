# Scientific workflow derived from Revision/RNAseq PCS met inhibitors/results/scripts/15_all_feature_strict_deg_figure_heatmaps_edgeR_exact.R
# Run from the repository root with METABOLIC_INPUT_DIR set to the external project data.
source("data_preprocessing/common/config.R")
suppressPackageStartupMessages({
  library(openxlsx)
  library(pheatmap)
})

analysis_root <- paste0(
  paste0(project_path("Revision/"), "/"),
  "RNAseq PCS met inhibitors/results/edgeR_batch_corrected_analysis"
)
qc_cache <- file.path(analysis_root, "cache", "01_qc_batch_correction.rds")
sensitivity_cache <- file.path(
  analysis_root,
  "cache",
  "05_method_threshold_sensitivity.rds"
)
output_dir <- Sys.getenv(
  "RNASEQ_FIGURE_OUTPUT_DIR",
  file.path(
    getwd(),
    "output",
    "heatmaps_v3_all_DEGs"
  )
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

qc <- readRDS(qc_cache)
sensitivity <- readRDS(sensitivity_cache)
expression_all <- qc$batch_corrected_logcpm
method_results <- sensitivity$all_method_results

fdr_cutoff <- 0.05
absolute_fc_cutoff <- 1.5
absolute_log2fc_cutoff <- log2(absolute_fc_cutoff)
z_limit <- 2.5

condition_colors <- c(
  "Control" = "#59636E",
  "2DG + OX" = "#C83E3A",
  "LDHA + OX" = "#8F4E93"
)
direction_colors <- c(
  "Downregulated" = "#2166AC",
  "Upregulated" = "#B2182B"
)

panel_profiles <- list(
  `90x90mm` = list(
    width_mm = 90,
    height_mm = 90,
    fontsize = 7.2,
    fontsize_col = 6.7,
    treeheight_row = 23
  ),
  `90x68mm` = list(
    width_mm = 90,
    height_mm = 68,
    fontsize = 6.5,
    fontsize_col = 6.1,
    treeheight_row = 19
  )
)

strict_features <- function(contrast_key) {
  selected <- method_results[
    as.character(method_results$method) == "edgeR_exact_ComBatSeq" &
      as.character(method_results$contrast) == contrast_key &
      is.finite(method_results$FDR) &
      is.finite(method_results$logFC) &
      method_results$FDR < fdr_cutoff &
      abs(method_results$logFC) > absolute_log2fc_cutoff,
    ,
    drop = FALSE
  ]
  selected <- selected[order(selected$FDR, -abs(selected$logFC)), , drop = FALSE]
  mapped <- !is.na(selected$SYMBOL) & nzchar(trimws(selected$SYMBOL))

  selected$Symbol <- ifelse(mapped, trimws(selected$SYMBOL), NA_character_)
  selected$Mapping_status <- factor(
    ifelse(mapped, "Mapped symbol", "No mapped symbol"),
    levels = c("Mapped symbol", "No mapped symbol")
  )
  selected$Direction <- factor(
    ifelse(selected$logFC < 0, "Downregulated", "Upregulated"),
    levels = c("Downregulated", "Upregulated")
  )
  selected$Fold_change <- 2^abs(selected$logFC)

  if (!nrow(selected)) stop("No strict feature-level DEGs for ", contrast_key)
  if (anyDuplicated(selected$gene_id)) {
    stop("Duplicate source feature IDs for ", contrast_key)
  }
  missing_expression <- setdiff(selected$gene_id, rownames(expression_all))
  if (length(missing_expression)) {
    stop(
      "Missing corrected expression for ",
      length(missing_expression),
      " strict features in ",
      contrast_key
    )
  }
  selected
}

sample_keys_for_heatmap <- function(treatment_code) {
  c(
    paste0("PCS", 1:5, "_C"),
    paste0("PCS", 1:5, "_", treatment_code)
  )
}

make_matrices <- function(features, treatment_code) {
  sample_keys <- sample_keys_for_heatmap(treatment_code)
  logcpm <- expression_all[features$gene_id, sample_keys, drop = FALSE]
  rownames(logcpm) <- features$gene_id

  row_means <- rowMeans(logcpm)
  row_sds <- apply(logcpm, 1, stats::sd)
  row_sds[!is.finite(row_sds) | row_sds == 0] <- 1
  z_score <- sweep(logcpm, 1, row_means, "-")
  z_score <- sweep(z_score, 1, row_sds, "/")

  correlation <- stats::cor(t(z_score), use = "pairwise.complete.obs")
  correlation[!is.finite(correlation)] <- 0
  diag(correlation) <- 1
  distance_matrix <- 1 - correlation
  distance_matrix[distance_matrix < 0] <- 0
  distance <- stats::as.dist(distance_matrix)
  row_tree <- stats::hclust(distance, method = "average")

  list(logCPM = logcpm, z_score = z_score, row_tree = row_tree)
}

clip_z_matrix <- function(x) {
  x[x < -z_limit] <- -z_limit
  x[x > z_limit] <- z_limit
  x
}

column_annotation <- function(treatment_code, treatment_label) {
  sample_keys <- sample_keys_for_heatmap(treatment_code)
  data.frame(
    Condition = factor(
      c(rep("Control", 5), rep(treatment_label, 5)),
      levels = c("Control", treatment_label)
    ),
    row.names = sample_keys
  )
}

row_annotation <- function(features) {
  data.frame(
    Direction = features$Direction,
    row.names = features$gene_id
  )
}

column_labels <- function(treatment_code) {
  c(
    paste0("C", 1:5),
    paste0(treatment_code, 1:5)
  )
}

build_heatmap <- function(
    matrices,
    features,
    treatment_code,
    treatment_label,
    title,
    profile) {
  pheatmap(
    clip_z_matrix(matrices$z_score),
    color = grDevices::colorRampPalette(
      c("#2166AC", "#F7F7F7", "#B2182B")
    )(101),
    breaks = seq(-z_limit, z_limit, length.out = 102),
    cluster_rows = matrices$row_tree,
    cluster_cols = FALSE,
    gaps_col = 5,
    annotation_row = row_annotation(features),
    annotation_col = column_annotation(treatment_code, treatment_label),
    annotation_colors = list(
      Condition = condition_colors[c("Control", treatment_label)],
      Direction = direction_colors
    ),
    annotation_names_row = FALSE,
    annotation_names_col = FALSE,
    show_rownames = FALSE,
    show_colnames = TRUE,
    labels_col = column_labels(treatment_code),
    angle_col = 45,
    border_color = NA,
    fontsize = profile$fontsize,
    fontsize_col = profile$fontsize_col,
    main = paste0(
      title,
      "\nFDR < 0.05; |FC| > 1.5",
      "\nPCS-corrected logCPM (row z-score)"
    ),
    legend_breaks = c(-2, 0, 2),
    legend_labels = c("-2", "0", "2"),
    treeheight_row = profile$treeheight_row,
    treeheight_col = 0,
    silent = TRUE
  )
}

save_pheatmap <- function(heatmap, stem, profile) {
  pdf_file <- file.path(output_dir, paste0(stem, ".pdf"))
  png_file <- file.path(output_dir, paste0(stem, ".png"))

  grDevices::cairo_pdf(
    filename = pdf_file,
    width = profile$width_mm / 25.4,
    height = profile$height_mm / 25.4,
    bg = "white"
  )
  grid::grid.newpage()
  grid::grid.draw(heatmap$gtable)
  grDevices::dev.off()

  grDevices::png(
    filename = png_file,
    width = profile$width_mm,
    height = profile$height_mm,
    units = "mm",
    res = 600,
    bg = "white",
    type = "cairo"
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
  features <- strict_features(contrast_key)
  matrices <- make_matrices(features, treatment_code)
  title <- paste0(
    treatment_label,
    " vs control: ",
    nrow(features),
    " DEGs"
  )

  files <- unlist(lapply(names(panel_profiles), function(profile_name) {
    profile <- panel_profiles[[profile_name]]
    heatmap <- build_heatmap(
      matrices = matrices,
      features = features,
      treatment_code = treatment_code,
      treatment_label = treatment_label,
      title = title,
      profile = profile
    )
    save_pheatmap(
      heatmap,
      paste0(file_prefix, "_", profile_name),
      profile
    )
  }))

  cluster_order <- data.frame(
    Cluster_rank = seq_len(nrow(features)),
    Feature_ID = rownames(matrices$z_score)[matrices$row_tree$order],
    stringsAsFactors = FALSE
  )
  features$Cluster_rank <- match(features$gene_id, cluster_order$Feature_ID)

  list(
    features = features,
    matrices = matrices,
    cluster_order = cluster_order,
    files = files
  )
}

result_2dg <- build_treatment_outputs(
  contrast_key = "2DG_vs_C",
  treatment_code = "2DG",
  treatment_label = "2DG + OX",
  file_prefix = "2DG_all_DEGs"
)
result_ldha <- build_treatment_outputs(
  contrast_key = "LDHA_vs_C",
  treatment_code = "LDHA",
  treatment_label = "LDHA + OX",
  file_prefix = "LDHA_all_DEGs"
)

expected_counts <- c(`2DG_vs_C` = 214L, `LDHA_vs_C` = 141L)
observed_counts <- c(
  `2DG_vs_C` = nrow(result_2dg$features),
  `LDHA_vs_C` = nrow(result_ldha$features)
)
if (!identical(observed_counts, expected_counts)) {
  stop(
    "Feature-level counts do not match the app: ",
    paste(names(observed_counts), observed_counts, collapse = ", ")
  )
}

feature_export <- function(result) {
  features <- result$features
  output <- data.frame(
    Cluster_rank = features$Cluster_rank,
    Feature_ID = features$gene_id,
    Symbol = features$Symbol,
    Mapping_status = as.character(features$Mapping_status),
    Direction = as.character(features$Direction),
    Contrast = features$contrast_label,
    log2FC = features$logFC,
    Fold_change = features$Fold_change,
    Statistic = features$statistic,
    P_value = features$PValue,
    FDR = features$FDR,
    Average_expression = features$average_expression,
    stringsAsFactors = FALSE
  )
  output[order(output$Cluster_rank), , drop = FALSE]
}

matrix_export <- function(matrix_values, cluster_order) {
  ordered <- matrix_values[cluster_order$Feature_ID, , drop = FALSE]
  data.frame(
    Cluster_rank = cluster_order$Cluster_rank,
    Feature_ID = rownames(ordered),
    ordered,
    check.names = FALSE,
    row.names = NULL
  )
}

features_2dg <- feature_export(result_2dg)
features_ldha <- feature_export(result_ldha)
csv_2dg <- file.path(
  output_dir,
  "2DG_all_feature_DEGs.csv"
)
csv_ldha <- file.path(
  output_dir,
  "LDHA_all_feature_DEGs.csv"
)
xlsx_file <- file.path(
  output_dir,
  "heatmap_values.xlsx"
)

write.csv(features_2dg, csv_2dg, row.names = FALSE, na = "")
write.csv(features_ldha, csv_ldha, row.names = FALSE, na = "")

notes_export <- data.frame(
  Item = c(
    "DEG method",
    "DEG threshold",
    "Heatmap rows",
    "Heatmap row labels",
    "Heatmap expression values",
    "Heatmap cell colors",
    "DEG direction annotation",
    "2DG + OX all feature-level DEGs",
    "2DG + OX mapped symbols",
    "2DG + OX features without mapped symbol",
    "2DG + OX all-feature direction",
    "LDHA + OX all feature-level DEGs",
    "LDHA + OX mapped symbols",
    "LDHA + OX features without mapped symbol",
    "LDHA + OX all-feature direction",
    "App up/down counters"
  ),
  Value = c(
    "edgeR exact test on ComBat-seq corrected counts",
    "FDR < 0.05 and absolute fold change > 1.5",
    "All feature-level DEGs, including features without a mapped symbol",
    "Hidden in every exported heatmap",
    "Row z-scores of ComBat-seq PCS-corrected logCPM",
    "Blue = lower relative expression; red = higher relative expression",
    "Blue = downregulated; red = upregulated",
    nrow(result_2dg$features),
    sum(result_2dg$features$Mapping_status == "Mapped symbol"),
    sum(result_2dg$features$Mapping_status == "No mapped symbol"),
    paste0(
      sum(result_2dg$features$Direction == "Upregulated"),
      " up; ",
      sum(result_2dg$features$Direction == "Downregulated"),
      " down"
    ),
    nrow(result_ldha$features),
    sum(result_ldha$features$Mapping_status == "Mapped symbol"),
    sum(result_ldha$features$Mapping_status == "No mapped symbol"),
    paste0(
      sum(result_ldha$features$Direction == "Upregulated"),
      " up; ",
      sum(result_ldha$features$Direction == "Downregulated"),
      " down"
    ),
    paste0(
      "The app reports up/down only among mapped symbols: ",
      "2DG 82/107 and LDHA 104/14"
    )
  ),
  stringsAsFactors = FALSE
)

write.xlsx(
  list(
    `README` = notes_export,
    `2DG_feature_DEGs` = features_2dg,
    `LDHA_feature_DEGs` = features_ldha,
    `2DG_logCPM` = matrix_export(
      result_2dg$matrices$logCPM,
      result_2dg$cluster_order
    ),
    `2DG_row_zscore` = matrix_export(
      result_2dg$matrices$z_score,
      result_2dg$cluster_order
    ),
    `LDHA_logCPM` = matrix_export(
      result_ldha$matrices$logCPM,
      result_ldha$cluster_order
    ),
    `LDHA_row_zscore` = matrix_export(
      result_ldha$matrices$z_score,
      result_ldha$cluster_order
    )
  ),
  xlsx_file,
  asTable = TRUE,
  overwrite = TRUE
)

manifest_file <- file.path(output_dir, "manifest.txt")
writeLines(
  c(
    "Figures: all-feature strict-DEG heatmaps for 2DG + OX and LDHA + OX",
    "Differential-expression method: edgeR exact test",
    "DEG definition: FDR < 0.05 and absolute fold change > 1.5",
    "Expression display: row z-scores of ComBat-seq PCS-corrected logCPM",
    "Heatmap cell colors: blue = lower relative expression; red = higher relative expression",
    "DEG direction annotation: blue = downregulated; red = upregulated",
    "Columns: ten individual samples; five controls followed by five treatment samples",
    "Rows: all called feature-level DEGs clustered by 1 - Pearson correlation using average linkage",
    "Row labels are hidden in all figure exports",
    "Columns are not clustered",
    paste0(
      "2DG + OX: 214 rows = 189 mapped symbols + 25 features without a mapped symbol; ",
      "98 up and 116 down across all features"
    ),
    paste0(
      "LDHA + OX: 141 rows = 118 mapped symbols + 23 features without a mapped symbol; ",
      "126 up and 15 down across all features"
    ),
    paste0(
      "The app up/down values remain symbol-only: ",
      "2DG 82 up/107 down; LDHA 104 up/14 down"
    ),
    "Figure presets: 90 x 90 mm and 90 x 68 mm; PDFs are vector and can be scaled in an A4 multipanel layout",
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
