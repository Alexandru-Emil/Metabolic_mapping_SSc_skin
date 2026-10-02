# Scientific workflow derived from Revision/RNAseq PCS met inhibitors/results/scripts/12_selected_deg_heatmaps_volcanoes_edgeR_exact.R
# Run from the repository root with METABOLIC_INPUT_DIR set to the external project data.
source("data_preprocessing/common/config.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
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
  file.path(getwd(), "output", "edgeR_exact_selected_DEG_heatmaps_volcanoes")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

report <- readRDS(report_file)
gene_results <- report$gene_results
expression_long <- report$expression
sample_labels <- report$sample_labels

heatmap_genes_2DG <- list(
  OXPHOS = c("NDUFV1", "COX3", "COX7A1"),
  Glycolysis = c("GALK1", "PFKL", "PFKM"),
  `Carbon catabolism` = c(
    "FAAH", "ECHDC2", "HIBADH", "LPIN3", "PLIN5", "LDHD"
  ),
  `ECM / fibrosis` = c(
    "CHI3L1", "COL9A2", "SFRP1", "FBLN2", "SERPINA3", "MMP7"
  )
)

volcano_labels_2DG <- c(
  "NDUFV1",
  "FAAH", "ECHDC2", "HIBADH", "LPIN3", "PLIN5", "LDHD",
  "CHI3L1", "COL9A2", "SFRP1", "FBLN2", "SERPINA3", "MMP7"
)

heatmap_genes_GNE <- list(
  OXPHOS = c("COX1", "ND5", "ND6", "COX3", "CYTB"),
  `Carbon catabolism` = c("ECHDC2", "IDUA", "DECR1"),
  `ECM / fibrosis` = c(
    "COL9A2", "MMP7", "SFRP1", "LOXL4", "SPARCL1"
  )
)

volcano_labels_GNE_strict <- c("COX1", "ND5", "ND6")
volcano_labels_GNE_contextual <- c("ECHDC2", "COL9A2", "COX3", "CYTB")

evidence_levels <- c(
  "Strict DEG",
  "Suggestive DEG",
  "GSEA leading-edge only"
)

pathway_colors <- c(
  "OXPHOS" = "#4477AA",
  "Glycolysis" = "#CCBB44",
  "Carbon catabolism" = "#228833",
  "ECM / fibrosis" = "#AA3377"
)
evidence_colors <- c(
  "Strict DEG" = "#B2182B",
  "Suggestive DEG" = "#D9A441",
  "GSEA leading-edge only" = "#8A949C"
)
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

flatten_gene_panel <- function(panel) {
  data.frame(
    Symbol = unlist(panel, use.names = FALSE),
    Pathway = rep(names(panel), lengths(panel)),
    Plot_order = seq_len(sum(lengths(panel))),
    stringsAsFactors = FALSE
  )
}

evidence_from_fdr <- function(fdr) {
  factor(
    ifelse(
      fdr < 0.05,
      "Strict DEG",
      ifelse(fdr < 0.10, "Suggestive DEG", "GSEA leading-edge only")
    ),
    levels = evidence_levels
  )
}

extract_panel_statistics <- function(contrast_key, panel) {
  annotation <- flatten_gene_panel(panel)
  contrast <- gene_results[
    gene_results$Contrast_key == contrast_key,
    ,
    drop = FALSE
  ]
  matches <- match(toupper(annotation$Symbol), toupper(contrast$Symbol))

  if (anyNA(matches)) {
    stop(
      "Missing requested symbols for ", contrast_key, ": ",
      paste(annotation$Symbol[is.na(matches)], collapse = ", ")
    )
  }

  selected <- contrast[matches, , drop = FALSE]
  if (anyDuplicated(toupper(selected$Symbol))) {
    stop("Duplicate symbols found in the selected panel for ", contrast_key, ".")
  }

  selected$Symbol <- annotation$Symbol
  selected$Pathway <- annotation$Pathway
  selected$Plot_order <- annotation$Plot_order
  selected$Evidence <- evidence_from_fdr(selected$FDR)
  selected[order(selected$Plot_order), , drop = FALSE]
}

panel_stats_2dg <- extract_panel_statistics("2DG_vs_C", heatmap_genes_2DG)
panel_stats_ldha <- extract_panel_statistics("LDHA_vs_C", heatmap_genes_GNE)

validate_label_set <- function(stats, symbols, label, require_strict = FALSE) {
  missing <- setdiff(symbols, stats$Symbol)
  if (length(missing)) {
    stop(label, " symbols missing from selected statistics: ", paste(missing, collapse = ", "))
  }
  if (require_strict) {
    non_strict <- stats$Symbol[
      stats$Symbol %in% symbols &
        !(stats$FDR < 0.05 & abs(stats$log2FC) >= log2(1.5))
    ]
    if (length(non_strict)) {
      stop(label, " contains non-strict labels: ", paste(non_strict, collapse = ", "))
    }
  }
}

validate_label_set(
  panel_stats_2dg,
  volcano_labels_2DG,
  "2DG volcano",
  require_strict = TRUE
)
validate_label_set(
  panel_stats_ldha,
  volcano_labels_GNE_strict,
  "LDHA strict volcano",
  require_strict = TRUE
)
validate_label_set(
  panel_stats_ldha,
  volcano_labels_GNE_contextual,
  "LDHA contextual volcano"
)

sample_keys_for_contrast <- function(treatment_code) {
  keys <- names(sample_labels)
  selected <- grepl("_C$", keys) | grepl(paste0("_", treatment_code, "$"), keys)
  keys[selected]
}

make_expression_matrices <- function(stats, treatment_code) {
  sample_keys <- sample_keys_for_contrast(treatment_code)
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
  row_index <- match(
    selected_expression$expression_key,
    stats$expression_key
  )
  col_index <- match(selected_expression$Sample_key, sample_keys)
  matrix_logcpm[cbind(row_index, col_index)] <- selected_expression$logCPM

  if (anyNA(matrix_logcpm)) {
    stop("Missing expression values in the ", treatment_code, " heatmap matrix.")
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

  list(logCPM = matrix_logcpm, z_score = matrix_z)
}

expression_2dg <- make_expression_matrices(panel_stats_2dg, "2DG")
expression_ldha <- make_expression_matrices(panel_stats_ldha, "LDHA")

make_column_annotation <- function(sample_keys, treatment_code, treatment_label) {
  condition <- ifelse(grepl("_C$", sample_keys), "Control", treatment_label)
  data.frame(
    Condition = factor(condition, levels = c("Control", treatment_label)),
    PCS = factor(sub("_.*$", "", sample_keys), levels = paste0("PCS", 1:5)),
    row.names = sample_keys,
    check.names = FALSE
  )
}

make_column_labels <- function(sample_keys, treatment_code) {
  condition_short <- ifelse(
    grepl("_C$", sample_keys),
    "C",
    ifelse(treatment_code == "2DG", "2DG", "LDHA")
  )
  paste(sub("_.*$", "", sample_keys), condition_short)
}

make_row_annotation <- function(stats) {
  data.frame(
    Pathway = factor(stats$Pathway, levels = unique(stats$Pathway)),
    Evidence = factor(stats$Evidence, levels = evidence_levels),
    row.names = stats$Symbol,
    check.names = FALSE
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

build_heatmap <- function(
    matrices,
    stats,
    panel,
    treatment_code,
    treatment_label,
    title) {
  sample_keys <- colnames(matrices$z_score)
  annotation_col <- make_column_annotation(
    sample_keys,
    treatment_code,
    treatment_label
  )
  annotation_row <- make_row_annotation(stats)
  heatmap_pathways <- unique(stats$Pathway)
  annotation_colors <- list(
    Condition = condition_colors[c("Control", treatment_label)],
    PCS = pcs_colors,
    Pathway = pathway_colors[heatmap_pathways],
    Evidence = evidence_colors
  )
  row_gaps <- head(cumsum(lengths(panel)), -1)
  column_gaps <- seq(2, ncol(matrices$z_score) - 2, by = 2)
  z_limit <- 2.5
  heatmap_matrix <- matrices$z_score
  heatmap_matrix[heatmap_matrix < -z_limit] <- -z_limit
  heatmap_matrix[heatmap_matrix > z_limit] <- z_limit

  pheatmap(
    heatmap_matrix,
    color = grDevices::colorRampPalette(
      c("#2166AC", "#F7F7F7", "#B2182B")
    )(101),
    breaks = seq(-z_limit, z_limit, length.out = 102),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    gaps_row = row_gaps,
    gaps_col = column_gaps,
    annotation_row = annotation_row,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    annotation_names_row = FALSE,
    annotation_names_col = TRUE,
    labels_col = make_column_labels(sample_keys, treatment_code),
    angle_col = 45,
    border_color = "#F1F3F4",
    fontsize = 8.2,
    fontsize_row = 8.2,
    fontsize_col = 7.3,
    main = paste0(
      title,
      "\nPCS-corrected logCPM (row z-score)"
    ),
    legend_breaks = c(-2, -1, 0, 1, 2),
    legend_labels = c("-2", "-1", "0", "1", "2"),
    treeheight_row = 0,
    treeheight_col = 0,
    silent = TRUE
  )
}

heatmap_2dg <- build_heatmap(
  matrices = expression_2dg,
  stats = panel_stats_2dg,
  panel = heatmap_genes_2DG,
  treatment_code = "2DG",
  treatment_label = "2DG + OX",
  title = "2DG + OX vs control"
)
heatmap_ldha <- build_heatmap(
  matrices = expression_ldha,
  stats = panel_stats_ldha,
  panel = heatmap_genes_GNE,
  treatment_code = "LDHA",
  treatment_label = "LDHA + OX",
  title = "LDHA + OX vs control"
)

heatmap_outputs <- c(
  save_pheatmap(
    heatmap_2dg,
    "2DG_OX_edgeR_exact_selected_gene_heatmap",
    width_mm = 165,
    height_mm = 165
  ),
  save_pheatmap(
    heatmap_ldha,
    "LDHA_OX_edgeR_exact_selected_gene_heatmap",
    width_mm = 165,
    height_mm = 145
  )
)

base_volcano_theme <- theme_bw(base_size = 9.2, base_family = "Arial") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#E4E7E9", linewidth = 0.28),
    panel.border = element_rect(color = "#4D555A", linewidth = 0.45),
    axis.text = element_text(color = "#202326"),
    axis.title = element_text(color = "#202326", size = 9),
    plot.title = element_text(face = "bold", size = 11, color = "#202326"),
    plot.subtitle = element_text(size = 8, color = "#4F5960"),
    plot.title.position = "plot",
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold", size = 7.5),
    legend.text = element_text(size = 7.2),
    plot.margin = margin(7, 10, 7, 7, unit = "pt")
  )

build_volcano <- function(
    contrast_key,
    title,
    strict_labels,
    contextual_labels = character()) {
  data <- gene_results[
    gene_results$Contrast_key == contrast_key &
      is.finite(gene_results$log2FC) &
      is.finite(gene_results$FDR),
    ,
    drop = FALSE
  ]
  fold_change_cutoff <- log2(1.5)
  data$minus_log10_FDR <- -log10(pmax(data$FDR, .Machine$double.xmin))
  data$Class <- "Not strict"
  data$Class[
    data$FDR < 0.05 & data$log2FC <= -fold_change_cutoff
  ] <- "Strict downregulated"
  data$Class[
    data$FDR < 0.05 & data$log2FC >= fold_change_cutoff
  ] <- "Strict upregulated"
  data$Class <- factor(
    data$Class,
    levels = c("Not strict", "Strict downregulated", "Strict upregulated")
  )
  data <- data[order(data$Class != "Not strict"), , drop = FALSE]

  strict <- data[data$Symbol %in% strict_labels, , drop = FALSE]
  strict <- strict[match(strict_labels, strict$Symbol), , drop = FALSE]
  contextual <- data[data$Symbol %in% contextual_labels, , drop = FALSE]
  contextual <- contextual[
    match(contextual_labels, contextual$Symbol),
    ,
    drop = FALSE
  ]
  if (anyNA(strict$Symbol) || anyNA(contextual$Symbol)) {
    stop("A requested volcano label was not found for ", contrast_key, ".")
  }
  strict_label_data <- strict
  strict_label_data$Volcano_label_type <- "Strict DEG label"
  contextual_label_data <- contextual
  contextual_label_data$Volcano_label_type <- rep(
    "Contextual label",
    nrow(contextual_label_data)
  )
  label_data <- rbind(strict_label_data, contextual_label_data)
  label_data$Volcano_label_type <- factor(
    label_data$Volcano_label_type,
    levels = c("Strict DEG label", "Contextual label")
  )

  plot <- ggplot(
    data,
    aes(x = log2FC, y = minus_log10_FDR, color = Class)
  ) +
    geom_point(size = 0.85, alpha = 0.60, stroke = 0) +
    geom_vline(
      xintercept = c(-fold_change_cutoff, fold_change_cutoff),
      linetype = "dashed",
      color = "#59636E",
      linewidth = 0.35
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      color = "#59636E",
      linewidth = 0.35
    ) +
    geom_point(
      data = strict,
      aes(x = log2FC, y = minus_log10_FDR),
      inherit.aes = FALSE,
      shape = 21,
      size = 2.5,
      fill = "#B2182B",
      color = "#202326",
      stroke = 0.45
    ) +
    scale_color_manual(
      name = "All-gene classification",
      values = c(
        "Not strict" = "#BFC5C9",
        "Strict downregulated" = "#B2182B",
        "Strict upregulated" = "#2166AC"
      ),
      drop = FALSE
    ) +
    labs(
      title = title,
      subtitle = paste0(
        "edgeR exact test; strict = FDR < 0.05 and |fold change| >= 1.5",
        if (length(contextual_labels)) {
          "; gold diamonds = contextual labels"
        } else {
          ""
        }
      ),
      x = "log2 fold change (treatment / control)",
      y = expression(-log[10] * "(edgeR FDR)")
    ) +
    guides(
      color = guide_legend(
        title.position = "top",
        nrow = 1,
        override.aes = list(size = 2.6, alpha = 1)
      )
    ) +
    coord_cartesian(clip = "off") +
    base_volcano_theme

  if (nrow(contextual)) {
    plot <- plot +
      geom_point(
        data = contextual,
        aes(x = log2FC, y = minus_log10_FDR),
        inherit.aes = FALSE,
        shape = 23,
        size = 2.7,
        fill = "#D9A441",
        color = "#202326",
        stroke = 0.48
      )
  }

  plot <- plot +
    ggrepel::geom_label_repel(
      data = label_data,
      aes(
        x = log2FC,
        y = minus_log10_FDR,
        label = Symbol,
        fill = Volcano_label_type
      ),
      inherit.aes = FALSE,
      seed = 20260726,
      size = 2.6,
      color = "#202326",
      label.size = 0.12,
      label.padding = grid::unit(0.10, "lines"),
      box.padding = 0.42,
      point.padding = 0.18,
      min.segment.length = 0,
      segment.color = "#6A7379",
      segment.size = 0.25,
      max.overlaps = Inf,
      max.time = 6,
      force = 11,
      show.legend = FALSE
    ) +
    scale_fill_manual(
      values = c(
        "Strict DEG label" = scales::alpha("white", 0.90),
        "Contextual label" = scales::alpha("#FFF1C2", 0.94)
      ),
      guide = "none"
    )

  list(plot = plot, all_genes = data, strict = strict, contextual = contextual)
}

volcano_2dg <- build_volcano(
  contrast_key = "2DG_vs_C",
  title = "2DG + OX vs control",
  strict_labels = volcano_labels_2DG
)
volcano_ldha <- build_volcano(
  contrast_key = "LDHA_vs_C",
  title = "LDHA + OX vs control",
  strict_labels = volcano_labels_GNE_strict,
  contextual_labels = volcano_labels_GNE_contextual
)

save_ggplot <- function(plot, stem, width_mm, height_mm) {
  pdf_file <- file.path(output_dir, paste0(stem, ".pdf"))
  png_file <- file.path(output_dir, paste0(stem, ".png"))

  ggsave(
    pdf_file,
    plot = plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    device = grDevices::cairo_pdf,
    bg = "white"
  )
  ggsave(
    png_file,
    plot = plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = 600,
    bg = "white"
  )
  c(pdf_file, png_file)
}

volcano_outputs <- c(
  save_ggplot(
    volcano_2dg$plot,
    "2DG_OX_edgeR_exact_labeled_volcano",
    width_mm = 150,
    height_mm = 132
  ),
  save_ggplot(
    volcano_ldha$plot,
    "LDHA_OX_edgeR_exact_labeled_volcano",
    width_mm = 145,
    height_mm = 125
  )
)

panel_export_columns <- c(
  "Symbol", "Pathway", "Evidence", "Contrast", "log2FC", "Fold_change",
  "Statistic", "P_value", "FDR", "Average_expression", "expression_key"
)
export_stats_2dg <- panel_stats_2dg[
  order(panel_stats_2dg$Plot_order),
  panel_export_columns
]
export_stats_ldha <- panel_stats_ldha[
  order(panel_stats_ldha$Plot_order),
  panel_export_columns
]

matrix_to_export <- function(matrix_values) {
  data.frame(
    Symbol = rownames(matrix_values),
    matrix_values,
    check.names = FALSE,
    row.names = NULL
  )
}

volcano_label_export <- rbind(
  data.frame(
    Treatment = "2DG + OX",
    Label_type = "Strict DEG label",
    volcano_2dg$strict[
      ,
      c("Symbol", "log2FC", "Fold_change", "P_value", "FDR"),
      drop = FALSE
    ],
    check.names = FALSE
  ),
  data.frame(
    Treatment = "LDHA + OX",
    Label_type = "Strict DEG label",
    volcano_ldha$strict[
      ,
      c("Symbol", "log2FC", "Fold_change", "P_value", "FDR"),
      drop = FALSE
    ],
    check.names = FALSE
  ),
  data.frame(
    Treatment = "LDHA + OX",
    Label_type = "Contextual label",
    volcano_ldha$contextual[
      ,
      c("Symbol", "log2FC", "Fold_change", "P_value", "FDR"),
      drop = FALSE
    ],
    check.names = FALSE
  )
)

csv_2dg <- file.path(output_dir, "2DG_OX_selected_gene_statistics.csv")
csv_ldha <- file.path(output_dir, "LDHA_OX_selected_gene_statistics.csv")
csv_labels <- file.path(output_dir, "edgeR_exact_volcano_label_statistics.csv")
xlsx_file <- file.path(output_dir, "edgeR_exact_selected_gene_plot_values.xlsx")

write.csv(export_stats_2dg, csv_2dg, row.names = FALSE, na = "")
write.csv(export_stats_ldha, csv_ldha, row.names = FALSE, na = "")
write.csv(volcano_label_export, csv_labels, row.names = FALSE, na = "")
write.xlsx(
  list(
    `2DG_gene_stats` = export_stats_2dg,
    `LDHA_gene_stats` = export_stats_ldha,
    `2DG_logCPM` = matrix_to_export(expression_2dg$logCPM),
    `2DG_row_zscore` = matrix_to_export(expression_2dg$z_score),
    `LDHA_logCPM` = matrix_to_export(expression_ldha$logCPM),
    `LDHA_row_zscore` = matrix_to_export(expression_ldha$z_score),
    `volcano_labels` = volcano_label_export
  ),
  xlsx_file,
  asTable = TRUE,
  overwrite = TRUE
)

manifest_file <- file.path(output_dir, "figure_manifest.txt")
writeLines(
  c(
    "Figures: selected-gene heatmaps and labeled volcano plots",
    "Differential-expression method: edgeR exact test",
    "Counts: ComBat-seq PCS-adjusted",
    "Heatmaps: row z-scores of corrected logCPM, Control and treatment only",
    "Heatmap sample order: paired Control/Treatment samples within PCS1-PCS5",
    "Heatmap rows: requested pathway panels in supplied order",
    "Evidence tiers: strict FDR < 0.05; suggestive 0.05 <= FDR < 0.10; otherwise GSEA leading-edge only",
    "Volcano axes: edgeR log2 fold change and -log10(FDR)",
    "Strict volcano background: FDR < 0.05 and absolute fold change >= 1.5",
    "LDHA contextual labels: gold diamonds, separate from strict DEG labels",
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    "2DG heatmap genes:",
    paste0(
      "- ",
      rep(names(heatmap_genes_2DG), lengths(heatmap_genes_2DG)),
      ": ",
      unlist(heatmap_genes_2DG, use.names = FALSE)
    ),
    "",
    "LDHA heatmap genes:",
    paste0(
      "- ",
      rep(names(heatmap_genes_GNE), lengths(heatmap_genes_GNE)),
      ": ",
      unlist(heatmap_genes_GNE, use.names = FALSE)
    )
  ),
  manifest_file
)

all_outputs <- c(
  heatmap_outputs,
  volcano_outputs,
  csv_2dg,
  csv_ldha,
  csv_labels,
  xlsx_file,
  manifest_file
)
cat("Created:\n", paste(all_outputs, collapse = "\n"), "\n")
