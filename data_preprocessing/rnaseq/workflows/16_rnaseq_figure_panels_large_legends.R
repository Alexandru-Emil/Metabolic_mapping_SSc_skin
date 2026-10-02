# Scientific workflow derived from Revision/RNAseq PCS met inhibitors/results/scripts/16_rnaseq_figure_panels_large_legends.R
# Run from the repository root with METABOLIC_INPUT_DIR set to the external project data.
source("data_preprocessing/common/config.R")
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(cowplot)
  library(ggplot2)
  library(svglite)
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
gsea_source_dir <- file.path(
  analysis_root,
  "publication_figures",
  "GSEA_edgeR_exact_merged_red_v2"
)
output_dir <- Sys.getenv(
  "RNASEQ_FIGURE_OUTPUT_DIR",
  file.path(getwd(), "output", "RNAseq_panels_v4_large_legends")
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

panel_width_mm <- 86
heatmap_height_mm <- 82
gsea_height_mm <- 108
composite_width_mm <- 178
composite_height_mm <- 205

condition_colors <- c(
  "Control" = "#59636E",
  "2DG + OX" = "#C83E3A",
  "LDHA + OX" = "#8F4E93"
)
direction_colors <- c(
  "Down" = "#2166AC",
  "Up" = "#B2182B"
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
  selected$Direction <- factor(
    ifelse(selected$logFC < 0, "Down", "Up"),
    levels = c("Down", "Up")
  )

  if (!nrow(selected)) stop("No strict feature-level DEGs for ", contrast_key)
  if (anyDuplicated(selected$gene_id)) stop("Duplicate feature IDs for ", contrast_key)
  if (!all(selected$gene_id %in% rownames(expression_all))) {
    stop("Corrected expression is missing for strict features in ", contrast_key)
  }
  selected
}

sample_keys <- function(treatment_code) {
  c(paste0("PCS", 1:5, "_C"), paste0("PCS", 1:5, "_", treatment_code))
}

heatmap_data <- function(contrast_key, treatment_code) {
  features <- strict_features(contrast_key)
  keys <- sample_keys(treatment_code)
  logcpm <- expression_all[features$gene_id, keys, drop = FALSE]
  rownames(logcpm) <- features$gene_id
  means <- rowMeans(logcpm)
  sds <- apply(logcpm, 1, stats::sd)
  sds[!is.finite(sds) | sds == 0] <- 1
  z <- sweep(sweep(logcpm, 1, means, "-"), 1, sds, "/")

  correlation <- stats::cor(t(z), use = "pairwise.complete.obs")
  correlation[!is.finite(correlation)] <- 0
  diag(correlation) <- 1
  distance_matrix <- 1 - correlation
  distance_matrix[distance_matrix < 0] <- 0
  row_tree <- stats::hclust(stats::as.dist(distance_matrix), method = "average")

  list(features = features, z = z, row_tree = row_tree)
}

legend_title_gp <- grid::gpar(fontsize = 9.5, fontface = "bold")
legend_label_gp <- grid::gpar(fontsize = 9)

make_heatmap <- function(
    contrast_key,
    treatment_code,
    treatment_label,
    title,
    expected_count) {
  data <- heatmap_data(contrast_key, treatment_code)
  if (nrow(data$z) != expected_count) {
    stop("Unexpected DEG count for ", contrast_key, ": ", nrow(data$z))
  }

  condition <- factor(
    c(rep("Control", 5), rep(treatment_label, 5)),
    levels = c("Control", treatment_label)
  )
  column_split <- factor(
    c(rep("Control", 5), rep("Treatment", 5)),
    levels = c("Control", "Treatment")
  )
  direction <- data$features$Direction
  names(direction) <- data$features$gene_id

  top_annotation <- HeatmapAnnotation(
    Condition = condition,
    col = list(Condition = condition_colors[c("Control", treatment_label)]),
    show_annotation_name = FALSE,
    simple_anno_size = grid::unit(3.2, "mm"),
    annotation_legend_param = list(
      Condition = list(
        title = "Condition",
        title_gp = legend_title_gp,
        labels_gp = legend_label_gp,
        grid_height = grid::unit(3.5, "mm"),
        grid_width = grid::unit(3.5, "mm"),
        gap = grid::unit(1.3, "mm")
      )
    )
  )
  left_annotation <- rowAnnotation(
    Direction = direction,
    col = list(Direction = direction_colors),
    show_annotation_name = FALSE,
    simple_anno_size = grid::unit(3.2, "mm"),
    annotation_legend_param = list(
      Direction = list(
        title = "DEG direction",
        title_gp = legend_title_gp,
        labels_gp = legend_label_gp,
        grid_height = grid::unit(3.5, "mm"),
        grid_width = grid::unit(3.5, "mm"),
        gap = grid::unit(1.3, "mm")
      )
    )
  )

  Heatmap(
    data$z,
    name = "Row z-score",
    col = circlize::colorRamp2(
      c(-z_limit, 0, z_limit),
      c("#2166AC", "#F7F7F7", "#B2182B")
    ),
    cluster_rows = data$row_tree,
    cluster_columns = FALSE,
    column_split = column_split,
    cluster_column_slices = FALSE,
    column_gap = grid::unit(2.2, "mm"),
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_labels = c(
      paste0("C", 1:5),
      paste0(treatment_code, 1:5)
    ),
    column_names_rot = 45,
    column_names_gp = grid::gpar(fontsize = 7.5),
    row_dend_width = grid::unit(10, "mm"),
    show_column_dend = FALSE,
    rect_gp = grid::gpar(col = NA),
    use_raster = FALSE,
    top_annotation = top_annotation,
    left_annotation = left_annotation,
    column_title = title,
    column_title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
    column_title_side = "top",
    heatmap_legend_param = list(
      title = "Row z-score",
      at = c(-2, 0, 2),
      labels = c("-2", "0", "2"),
      title_gp = legend_title_gp,
      labels_gp = legend_label_gp,
      grid_width = grid::unit(4.2, "mm"),
      legend_height = grid::unit(26, "mm")
    )
  )
}

capture_heatmap <- function(heatmap, width_mm, height_mm) {
  grid::grid.grabExpr(
    draw(
      heatmap,
      newpage = FALSE,
      merge_legends = TRUE,
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      legend_gap = grid::unit(2.5, "mm"),
      padding = grid::unit(c(2.5, 2.5, 2.5, 2.5), "mm")
    ),
    width = width_mm / 25.4,
    height = height_mm / 25.4
  )
}

save_grob <- function(grob, stem, width_mm, height_mm) {
  pdf_file <- file.path(output_dir, paste0(stem, ".pdf"))
  svg_file <- file.path(output_dir, paste0(stem, ".svg"))
  png_file <- file.path(output_dir, paste0(stem, ".png"))

  grDevices::cairo_pdf(
    pdf_file,
    width = width_mm / 25.4,
    height = height_mm / 25.4,
    bg = "white"
  )
  grid::grid.newpage()
  grid::grid.draw(grob)
  grDevices::dev.off()

  svglite::svglite(
    svg_file,
    width = width_mm / 25.4,
    height = height_mm / 25.4,
    bg = "white"
  )
  grid::grid.newpage()
  grid::grid.draw(grob)
  grDevices::dev.off()

  grDevices::png(
    png_file,
    width = width_mm,
    height = height_mm,
    units = "mm",
    res = 600,
    bg = "white",
    type = "cairo"
  )
  grid::grid.newpage()
  grid::grid.draw(grob)
  grDevices::dev.off()

  c(pdf_file, svg_file, png_file)
}

heatmap_2dg <- make_heatmap(
  contrast_key = "2DG_vs_C",
  treatment_code = "2DG",
  treatment_label = "2DG + OX",
  title = "2DG + OX vs control (214 DEGs)",
  expected_count = 214L
)
heatmap_ldha <- make_heatmap(
  contrast_key = "LDHA_vs_C",
  treatment_code = "LDHA",
  treatment_label = "LDHA + OX",
  title = "LDHA + OX vs control (141 DEGs)",
  expected_count = 141L
)

heatmap_grob_2dg <- capture_heatmap(
  heatmap_2dg,
  panel_width_mm,
  heatmap_height_mm
)
heatmap_grob_ldha <- capture_heatmap(
  heatmap_ldha,
  panel_width_mm,
  heatmap_height_mm
)

read_gsea_panel_data <- function(filename, treatment_order) {
  data <- read.csv(
    file.path(gsea_source_dir, filename),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  expected <- c("Category", "Term", "NES", "Set_size", "FDR")
  if (!all(expected %in% names(data))) {
    stop("Missing required GSEA columns in ", filename)
  }

  category_levels <- c(
    "OXPHOS",
    "Glycolysis / carbon catabolism",
    "ECM / fibrosis"
  )
  data$Category <- factor(data$Category, levels = category_levels)
  wrapped <- vapply(
    treatment_order,
    function(term) paste(strwrap(term, width = 29), collapse = "\n"),
    character(1)
  )
  names(wrapped) <- treatment_order
  data$Term_label <- factor(
    unname(wrapped[data$Term]),
    levels = rev(unname(wrapped[treatment_order]))
  )
  data$minus_log10_FDR <- -log10(pmax(data$FDR, .Machine$double.xmin))
  data$FDR_status <- ifelse(data$FDR <= 0.25, "FDR <= 0.25", "FDR > 0.25")
  data
}

terms_2dg <- c(
  "respiratory chain complex",
  "mitochondrial ATP synthesis coupled electron transport",
  "oxidative phosphorylation",
  "aerobic respiration",
  "hexose catabolic process",
  "small molecule catabolic process",
  "organic acid catabolic process",
  "monocarboxylic acid catabolic process",
  "extracellular matrix",
  "collagen-containing extracellular matrix",
  "extracellular matrix structural constituent"
)
terms_ldha <- c(
  "aerobic electron transport chain",
  "mitochondrial ATP synthesis coupled electron transport",
  "respiratory chain complex",
  "oxidative phosphorylation",
  "small molecule catabolic process",
  "organic acid catabolic process",
  "monocarboxylic acid catabolic process",
  "collagen fibril organization",
  "extracellular matrix structural constituent",
  "collagen-containing extracellular matrix",
  "extracellular matrix"
)

gsea_2dg <- read_gsea_panel_data(
  "2DG_OX_edgeR_exact_selected_GSEA_values.csv",
  terms_2dg
)
gsea_ldha <- read_gsea_panel_data(
  "LDHA_OX_edgeR_exact_selected_GSEA_values.csv",
  terms_ldha
)

all_nes <- c(gsea_2dg$NES, gsea_ldha$NES)
x_limits <- c(
  floor(min(c(0, all_nes), na.rm = TRUE) * 10) / 10 - 0.08,
  ceiling(max(c(0, all_nes), na.rm = TRUE) * 10) / 10 + 0.08
)
x_breaks <- c(-2.5, -1.5, 0)
size_limits <- range(c(gsea_2dg$Set_size, gsea_ldha$Set_size), na.rm = TRUE)
size_breaks <- c(100, 250, 400)

make_gsea_plot <- function(data, title) {
  significant <- data[data$FDR <= 0.25, , drop = FALSE]
  nonsignificant <- data[data$FDR > 0.25, , drop = FALSE]
  significant_scale <- significant$minus_log10_FDR
  if (!length(significant_scale)) stop("No GSEA terms at FDR <= 0.25 for ", title)
  fdr_limits <- range(significant_scale, na.rm = TRUE)
  if (diff(fdr_limits) == 0) fdr_limits <- fdr_limits + c(-0.05, 0.05)
  fdr_breaks <- seq(fdr_limits[1], fdr_limits[2], length.out = 3)

  plot <- ggplot(data, aes(x = NES, y = Term_label, size = Set_size)) +
    geom_point(
      data = significant,
      aes(fill = minus_log10_FDR),
      shape = 21,
      color = "#202326",
      stroke = 0.45,
      alpha = 0.98
    )
  if (nrow(nonsignificant)) {
    plot <- plot +
      geom_point(
        data = nonsignificant,
        aes(color = FDR_status),
        shape = 21,
        fill = "#B8BDC3",
        stroke = 0.45,
        alpha = 0.98
      )
  }

  guides_list <- list(
    fill = guide_colorbar(
      order = 1,
      title.position = "top",
      title.theme = element_text(size = 9, face = "bold"),
      label.theme = element_text(size = 8.3),
      barwidth = grid::unit(25, "mm"),
      barheight = grid::unit(3.2, "mm")
    ),
    size = guide_legend(
      order = 2,
      title.position = "top",
      title.theme = element_text(size = 9, face = "bold"),
      label.theme = element_text(size = 8.3),
      nrow = 1,
      override.aes = list(shape = 21, fill = "#F4B3B0", color = "#202326")
    )
  )
  if (nrow(nonsignificant)) {
    guides_list$colour <- guide_legend(
      order = 3,
      title.position = "top",
      title.theme = element_text(size = 9, face = "bold"),
      label.theme = element_text(size = 8.3),
      nrow = 1,
      override.aes = list(
        shape = 21,
        fill = "#B8BDC3",
        colour = "#202326",
        size = 3.8
      )
    )
  }

  plot <- plot +
    geom_vline(
      xintercept = 0,
      color = "#4D555A",
      linewidth = 0.4,
      linetype = "dashed"
    ) +
    scale_x_continuous(
      name = "NES",
      limits = x_limits,
      breaks = x_breaks,
      expand = expansion(mult = c(0.03, 0.03))
    ) +
    scale_size_continuous(
      name = "GO gene-set size",
      range = c(2.5, 5.8),
      limits = size_limits,
      breaks = size_breaks
    ) +
    scale_fill_gradient(
      name = expression(-log[10] * "(FDR)"),
      low = "#FADBD8",
      high = "#B2182B",
      limits = fdr_limits,
      breaks = fdr_breaks,
      labels = scales::label_number(accuracy = 0.1)
    )

  if (nrow(nonsignificant)) {
    plot <- plot +
      scale_color_manual(
        name = "FDR status",
        values = c("FDR > 0.25" = "#202326"),
        drop = TRUE
      )
  }

  plot +
    facet_grid(
      rows = vars(Category),
      scales = "free_y",
      space = "free_y",
      switch = "y",
      drop = TRUE,
      labeller = as_labeller(
        c(
          "OXPHOS" = "OXPHOS",
          "Glycolysis / carbon catabolism" = "Glycolysis /\ncarbon\ncatabolism",
          "ECM / fibrosis" = "ECM /\nfibrosis"
        )
      )
    ) +
    labs(title = title, y = NULL) +
    do.call(guides, guides_list) +
    coord_cartesian(clip = "off") +
    theme_bw(base_size = 9, base_family = "sans") +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "#D9DEE2", linewidth = 0.3),
      panel.border = element_rect(color = "#4D555A", linewidth = 0.5),
      strip.background = element_rect(
        fill = "#EEF1F3",
        color = "#9FAAB1",
        linewidth = 0.45
      ),
      strip.placement = "outside",
      strip.text.y.left = element_text(
        face = "bold",
        color = "#202326",
        angle = 0,
        size = 6.8,
        lineheight = 0.95
      ),
      axis.text.x = element_text(size = 7.5, color = "#202326"),
      axis.text.y = element_text(size = 6.6, color = "#202326", lineheight = 0.94),
      axis.title.x = element_text(size = 8.5, color = "#202326"),
      plot.title = element_text(
        face = "bold",
        size = 9.2,
        color = "#202326",
        hjust = 0,
        margin = margin(l = 14, unit = "pt")
      ),
      plot.title.position = "plot",
      legend.position = "bottom",
      legend.box = "vertical",
      legend.direction = "horizontal",
      legend.box.just = "center",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8.3),
      legend.key.height = grid::unit(4.2, "mm"),
      legend.key.width = grid::unit(5.2, "mm"),
      legend.spacing.y = grid::unit(1.2, "mm"),
      legend.spacing.x = grid::unit(1.8, "mm"),
      legend.margin = margin(t = 1.5, unit = "pt"),
      plot.margin = margin(4, 4, 3, 4, unit = "pt")
    )
}

gsea_plot_2dg <- make_gsea_plot(gsea_2dg, "GSEA: 2DG + OX vs control")
gsea_plot_ldha <- make_gsea_plot(gsea_ldha, "GSEA: LDHA + OX vs control")

save_ggplot <- function(plot, stem, width_mm, height_mm) {
  pdf_file <- file.path(output_dir, paste0(stem, ".pdf"))
  svg_file <- file.path(output_dir, paste0(stem, ".svg"))
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
    svg_file,
    plot = plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    device = svglite::svglite,
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
  c(pdf_file, svg_file, png_file)
}

heatmap_plot_2dg <- cowplot::ggdraw() + cowplot::draw_grob(heatmap_grob_2dg)
heatmap_plot_ldha <- cowplot::ggdraw() + cowplot::draw_grob(heatmap_grob_ldha)
top_row <- cowplot::plot_grid(
  heatmap_plot_2dg,
  heatmap_plot_ldha,
  ncol = 2,
  labels = c("B", "C"),
  label_size = 11,
  label_fontface = "plain",
  label_x = 0,
  label_y = 1,
  hjust = -0.15,
  vjust = 1.1
)
bottom_row <- cowplot::plot_grid(
  gsea_plot_2dg,
  gsea_plot_ldha,
  ncol = 2,
  labels = c("D", "E"),
  label_size = 10,
  label_fontface = "plain",
  label_x = 0,
  label_y = 1,
  hjust = 0,
  vjust = 1.1
)
composite <- cowplot::plot_grid(
  top_row,
  bottom_row,
  ncol = 1,
  rel_heights = c(heatmap_height_mm, gsea_height_mm)
)

outputs <- c(
  save_grob(
    heatmap_grob_2dg,
    "B_2DG_heatmap_large_legends",
    panel_width_mm,
    heatmap_height_mm
  ),
  save_grob(
    heatmap_grob_ldha,
    "C_LDHA_heatmap_large_legends",
    panel_width_mm,
    heatmap_height_mm
  ),
  save_ggplot(
    gsea_plot_2dg,
    "D_2DG_GSEA_large_legends",
    panel_width_mm,
    gsea_height_mm
  ),
  save_ggplot(
    gsea_plot_ldha,
    "E_LDHA_GSEA_large_legends",
    panel_width_mm,
    gsea_height_mm
  ),
  save_ggplot(
    composite,
    "BCDE_RNAseq_block_large_legends",
    composite_width_mm,
    composite_height_mm
  )
)

legend_file <- file.path(output_dir, "figure_legend.txt")
writeLines(
  c(
    paste0(
      "B-C, Differential-expression heatmaps from edgeR exact tests on ComBat-seq ",
      "PCS-corrected counts. Rows comprise all feature-level DEGs at FDR < 0.05 ",
      "and absolute fold change > 1.5, including features without mapped symbols ",
      "(214 for 2DG + OX and 141 for LDHA + OX). Columns are individual control ",
      "and treatment samples. Values are PCS-corrected logCPM transformed to ",
      "row z-scores; rows were clustered by 1 - Pearson correlation with average linkage. ",
      "Red denotes higher relative expression or upregulation, and blue denotes lower ",
      "relative expression or downregulation."
    ),
    paste0(
      "D-E, GO GSEA dot plots from the same edgeR exact analysis. Dot position is the ",
      "normalized enrichment score (NES), size is the GO gene-set size, and red intensity ",
      "is -log10(FDR). The red color range is scaled separately within each treatment. ",
      "Terms with FDR > 0.25, if present, are gray."
    )
  ),
  legend_file
)

manifest_file <- file.path(output_dir, "manifest.txt")
writeLines(
  c(
    "RNA-seq replacement panels B-E with enlarged final-size legends",
    "Analysis and term selection are unchanged from the existing edgeR exact outputs",
    "Individual panel width: 86 mm",
    "Heatmap panel height: 82 mm",
    "GSEA panel height: 108 mm",
    "Combined B-E block: 178 x 205 mm",
    "Heatmap legend text: 9 pt; title: 11 pt; sample labels: 7.5 pt",
    "GSEA legend titles: 9 pt; legend labels: 8.3 pt; title: 9.2 pt",
    "Exports: vector PDF, vector SVG, and 600 dpi PNG",
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
  ),
  manifest_file
)

cat("Created:\n", paste(c(outputs, legend_file, manifest_file), collapse = "\n"), "\n")
