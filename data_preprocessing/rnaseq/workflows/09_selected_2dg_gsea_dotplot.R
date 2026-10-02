# Scientific workflow derived from Revision/RNAseq PCS met inhibitors/results/scripts/09_selected_2dg_gsea_dotplot.R
# Run from the repository root with METABOLIC_INPUT_DIR set to the external project data.
source("data_preprocessing/common/config.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(openxlsx)
})

report_root <- project_path("Revision/RNAseq PCS met inhibitors/results/edgeR_batch_corrected_analysis/report")
output_dir <- Sys.getenv(
  "RNASEQ_FIGURE_OUTPUT_DIR",
  file.path(getwd(), "output", "2DG_GSEA_selected_terms")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

terms_to_plot <- c(
  "respiratory chain complex",
  "mitochondrial ATP synthesis coupled electron transport",
  "oxidative phosphorylation",
  "aerobic respiration",
  "mitochondrial inner membrane",
  "hexose catabolic process",
  "glycolytic process through glucose-6-phosphate",
  "extracellular matrix",
  "collagen-containing extracellular matrix"
)

term_annotation <- data.frame(
  Term = terms_to_plot,
  Category = c(rep("OXPHOS", 5), rep("Glycolysis", 2), rep("ECM", 2)),
  Plot_order = seq_along(terms_to_plot),
  stringsAsFactors = FALSE
)

report_specs <- list(
  "edgeR exact" = "edgeR_exact_ComBatSeq_FDR05_FC15_threshold_GO_GSEA_v5/report_data.rds",
  "DESeq2" = "DESeq2_ComBatSeq_FDR05_FC15_threshold_GO_GSEA_v5/report_data.rds"
)

extract_terms <- function(method, report_path) {
  report <- readRDS(file.path(report_root, report_path))
  data <- report$gsea_results
  data <- data[data$Contrast_key == "2DG_vs_C", , drop = FALSE]
  selected <- data[tolower(data$Term) %in% tolower(terms_to_plot), , drop = FALSE]

  matches <- match(tolower(terms_to_plot), tolower(selected$Term))
  if (anyNA(matches)) {
    stop("Missing terms for ", method, ": ", paste(terms_to_plot[is.na(matches)], collapse = ", "))
  }
  selected <- selected[matches, , drop = FALSE]
  if (anyDuplicated(tolower(selected$Term))) {
    stop("Duplicate requested terms found for ", method, ".")
  }

  selected$Method <- method
  selected <- merge(selected, term_annotation, by = "Term", sort = FALSE)
  selected[order(selected$Plot_order), , drop = FALSE]
}

plot_data <- do.call(
  rbind,
  Map(extract_terms, names(report_specs), unname(report_specs))
)
rownames(plot_data) <- NULL

plot_data$Method <- factor(plot_data$Method, levels = names(report_specs))
plot_data$Category <- factor(plot_data$Category, levels = c("OXPHOS", "Glycolysis", "ECM"))
term_labels <- setNames(
  vapply(
    terms_to_plot,
    function(term) paste(strwrap(term, width = 38), collapse = "\n"),
    character(1)
  ),
  terms_to_plot
)
plot_data$Term_label <- factor(
  unname(term_labels[plot_data$Term]),
  levels = rev(unname(term_labels[terms_to_plot]))
)
plot_data$minus_log10_FDR <- -log10(pmax(plot_data$FDR, .Machine$double.xmin))
plot_data$FDR_status <- ifelse(plot_data$FDR < 0.05, "FDR < 0.05", "FDR >= 0.05")

common_theme <- function(base_size = 9) {
  theme_bw(base_size = base_size, base_family = "Arial") +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "#D9DEE2", linewidth = 0.28),
      panel.border = element_rect(color = "#4D555A", linewidth = 0.45),
      strip.background = element_rect(fill = "#EEF1F3", color = "#9FAAB1", linewidth = 0.4),
      strip.text.y.left = element_text(face = "bold", color = "#202326", angle = 0, size = 7.5),
      axis.text = element_text(color = "#202326"),
      axis.text.y = element_text(size = 7.8),
      axis.title = element_text(color = "#202326", size = 8.5),
      plot.title = element_text(face = "bold", size = 10.5, color = "#202326"),
      plot.subtitle = element_text(size = 8, color = "#4F5960"),
      plot.title.position = "plot",
      legend.position = "bottom",
      legend.box = "vertical",
      legend.direction = "horizontal",
      legend.title = element_text(size = 7.5, face = "bold"),
      legend.text = element_text(size = 7),
      legend.key.height = grid::unit(3.5, "mm"),
      legend.key.width = grid::unit(5, "mm"),
      legend.spacing.x = grid::unit(2, "mm"),
      legend.margin = margin(t = -4, unit = "pt"),
      plot.margin = margin(7, 8, 6, 7, unit = "pt")
    )
}

scale_and_labels <- function(subtitle) {
  list(
    geom_vline(xintercept = 0, color = "#4D555A", linewidth = 0.35, linetype = "dashed"),
    scale_x_continuous(
      name = "NES",
      limits = c(floor(min(plot_data$NES) * 10) / 10 - 0.05, 0.05),
      breaks = seq(-3, 0, by = 0.5),
      expand = expansion(mult = c(0.04, 0.01))
    ),
    scale_size_continuous(
      name = "Gene-set size",
      range = c(2.6, 6),
      breaks = c(100, 300, 450)
    ),
    scale_fill_gradient(
      name = expression(-log[10] * "(FDR)"),
      low = "#DCE5EA",
      high = "#D9822B"
    ),
    facet_grid(
      rows = vars(Category),
      scales = "free_y",
      space = "free_y",
      switch = "y"
    ),
    labs(
      title = "2DG + OX vs control",
      subtitle = subtitle,
      y = NULL
    ),
    guides(
      fill = guide_colorbar(
        order = 1,
        title.position = "top",
        barwidth = grid::unit(25, "mm"),
        barheight = grid::unit(2.8, "mm")
      ),
      size = guide_legend(order = 2, title.position = "top", nrow = 1),
      shape = guide_legend(order = 3, title.position = "top", nrow = 1)
    ),
    coord_cartesian(clip = "off"),
    common_theme()
  )
}

comparison_plot <- ggplot(
  plot_data,
  aes(
    x = NES,
    y = Term_label,
    size = Set_size,
    fill = minus_log10_FDR,
    shape = Method,
    group = Method
  )
) +
  geom_point(
    position = position_dodge(width = 0.46),
    color = "#202326",
    stroke = 0.45,
    alpha = 0.96
  ) +
  scale_shape_manual(
    name = "DE method",
    values = c("edgeR exact" = 21, "DESeq2" = 24)
  ) +
  scale_and_labels("Selected GO GSEA terms; edgeR exact and DESeq2")

method_plot <- function(method) {
  method_data <- plot_data[plot_data$Method == method, , drop = FALSE]
  ggplot(
    method_data,
    aes(x = NES, y = Term_label, size = Set_size, fill = minus_log10_FDR)
  ) +
    geom_point(shape = 21, color = "#202326", stroke = 0.45, alpha = 0.96) +
    scale_and_labels(paste("Selected GO GSEA terms;", method))
}

save_plot <- function(plot, stem, width_mm, height_mm) {
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

outputs <- c(
  save_plot(
    comparison_plot,
    "2DG_OX_GSEA_selected_terms_edgeR_vs_DESeq2_dotplot",
    width_mm = 150,
    height_mm = 145
  ),
  save_plot(
    method_plot("edgeR exact"),
    "2DG_OX_GSEA_selected_terms_edgeR_exact_dotplot",
    width_mm = 125,
    height_mm = 145
  ),
  save_plot(
    method_plot("DESeq2"),
    "2DG_OX_GSEA_selected_terms_DESeq2_dotplot",
    width_mm = 125,
    height_mm = 145
  )
)

export_columns <- c(
  "Method", "Category", "Plot_order", "GO_ID", "Term", "Ontology",
  "Direction", "NES", "Enrichment_score", "Set_size", "Core_gene_count",
  "P_value", "FDR", "q_value", "FDR_status", "Core_genes"
)
export_data <- plot_data[order(plot_data$Plot_order, plot_data$Method), export_columns]

csv_file <- file.path(output_dir, "2DG_OX_GSEA_selected_terms_values.csv")
xlsx_file <- file.path(output_dir, "2DG_OX_GSEA_selected_terms_values.xlsx")
write.csv(export_data, csv_file, row.names = FALSE, na = "")
write.xlsx(
  list(
    comparison = export_data,
    edgeR_exact = export_data[export_data$Method == "edgeR exact", , drop = FALSE],
    DESeq2 = export_data[export_data$Method == "DESeq2", , drop = FALSE]
  ),
  xlsx_file,
  asTable = TRUE,
  overwrite = TRUE
)

manifest_file <- file.path(output_dir, "figure_manifest.txt")
writeLines(
  c(
    "Figure: 2DG + OX versus control selected GO GSEA terms",
    "Ranking: sign(log2FC) x -log10(raw P value)",
    "Counts: ComBat-seq PCS-corrected",
    "X axis: normalized enrichment score (NES)",
    "Dot size: GO gene-set size",
    "Dot fill: -log10(FDR)",
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    "Requested terms:",
    paste0("- ", term_annotation$Category, ": ", term_annotation$Term)
  ),
  manifest_file
)

cat("Created:\n", paste(c(outputs, csv_file, xlsx_file, manifest_file), collapse = "\n"), "\n")
