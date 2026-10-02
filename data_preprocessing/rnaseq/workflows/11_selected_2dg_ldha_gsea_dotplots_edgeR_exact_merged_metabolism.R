# Scientific workflow derived from Revision/RNAseq PCS met inhibitors/results/scripts/11_selected_2dg_ldha_gsea_dotplots_edgeR_exact_merged_metabolism.R
# Run from the repository root with METABOLIC_INPUT_DIR set to the external project data.
source("data_preprocessing/common/config.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(openxlsx)
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
    "2DG_LDHA_GSEA_selected_terms_edgeR_exact_merged_metabolism"
  )
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c(
  "OXPHOS",
  "Glycolysis / carbon catabolism",
  "ECM / fibrosis"
)

terms_2DG_IACS <- c(
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

groups_2DG_IACS <- c(
  rep("OXPHOS", 4),
  rep("Glycolysis / carbon catabolism", 4),
  rep("ECM / fibrosis", 3)
)

terms_GNE_IACS <- c(
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

groups_GNE_IACS <- c(
  rep("OXPHOS", 4),
  rep("Glycolysis / carbon catabolism", 3),
  rep("ECM / fibrosis", 4)
)

if (length(terms_2DG_IACS) != length(groups_2DG_IACS)) {
  stop("The 2DG term and group vectors have different lengths.")
}
if (length(terms_GNE_IACS) != length(groups_GNE_IACS)) {
  stop("The LDHA term and group vectors have different lengths.")
}

term_groups <- rbind(
  data.frame(
    Term = terms_2DG_IACS,
    Category = groups_2DG_IACS,
    Source_list = "2DG + OX",
    Source_order = seq_along(terms_2DG_IACS),
    stringsAsFactors = FALSE
  ),
  data.frame(
    Term = terms_GNE_IACS,
    Category = groups_GNE_IACS,
    Source_list = "LDHA + OX",
    Source_order = seq_along(terms_GNE_IACS),
    stringsAsFactors = FALSE
  )
)

category_by_term <- unique(term_groups[c("Term", "Category")])
conflicted_terms <- names(which(table(category_by_term$Term) > 1))
if (length(conflicted_terms)) {
  stop(
    "Terms assigned to more than one group: ",
    paste(conflicted_terms, collapse = ", ")
  )
}

union_terms <- unique(c(terms_2DG_IACS, terms_GNE_IACS))
union_annotation <- category_by_term[
  match(union_terms, category_by_term$Term),
  ,
  drop = FALSE
]
union_annotation$Union_order <- seq_len(nrow(union_annotation))

report <- readRDS(report_file)
gsea <- report$gsea_results

extract_exact_terms <- function(contrast_key, treatment, terms, groups) {
  contrast_data <- gsea[gsea$Contrast_key == contrast_key, , drop = FALSE]
  matches <- match(tolower(terms), tolower(contrast_data$Term))

  if (anyNA(matches)) {
    stop(
      "Missing terms for ", treatment, ": ",
      paste(terms[is.na(matches)], collapse = ", ")
    )
  }

  selected <- contrast_data[matches, , drop = FALSE]
  if (anyDuplicated(tolower(selected$Term))) {
    stop("Duplicate requested terms found for ", treatment, ".")
  }

  selected$Treatment <- treatment
  selected$Category <- groups
  selected$Treatment_order <- seq_along(terms)
  selected
}

extract_union_terms <- function(contrast_key, treatment) {
  selected <- extract_exact_terms(
    contrast_key = contrast_key,
    treatment = treatment,
    terms = union_terms,
    groups = union_annotation$Category
  )
  selected$Union_order <- union_annotation$Union_order
  selected
}

data_2dg <- extract_exact_terms(
  contrast_key = "2DG_vs_C",
  treatment = "2DG + OX",
  terms = terms_2DG_IACS,
  groups = groups_2DG_IACS
)
data_ldha <- extract_exact_terms(
  contrast_key = "LDHA_vs_C",
  treatment = "LDHA + OX",
  terms = terms_GNE_IACS,
  groups = groups_GNE_IACS
)
combined_data <- rbind(
  extract_union_terms("2DG_vs_C", "2DG + OX"),
  extract_union_terms("LDHA_vs_C", "LDHA + OX")
)
rownames(combined_data) <- NULL

prepare_plot_data <- function(data, term_order) {
  data$Treatment <- factor(
    data$Treatment,
    levels = c("2DG + OX", "LDHA + OX")
  )
  data$Category <- factor(data$Category, levels = group_levels)

  term_labels <- setNames(
    vapply(
      term_order,
      function(term) paste(strwrap(term, width = 39), collapse = "\n"),
      character(1)
    ),
    term_order
  )
  data$Term_label <- factor(
    unname(term_labels[data$Term]),
    levels = rev(unname(term_labels[term_order]))
  )
  data$minus_log10_FDR <- -log10(pmax(data$FDR, .Machine$double.xmin))
  data$FDR_status <- ifelse(
    data$FDR <= 0.25,
    "FDR <= 0.25",
    "FDR > 0.25"
  )
  data
}

plot_2dg <- prepare_plot_data(data_2dg, terms_2DG_IACS)
plot_ldha <- prepare_plot_data(data_ldha, terms_GNE_IACS)
plot_combined <- prepare_plot_data(combined_data, union_terms)

all_nes <- plot_combined$NES
x_limits <- c(
  floor(min(c(0, all_nes), na.rm = TRUE) * 10) / 10 - 0.08,
  ceiling(max(c(0, all_nes), na.rm = TRUE) * 10) / 10 + 0.08
)
x_breaks <- pretty(x_limits, n = 6)
size_limits <- range(plot_combined$Set_size, na.rm = TRUE)
size_breaks <- scales::breaks_pretty(n = 3)(size_limits)

common_theme <- function(base_size = 9) {
  theme_bw(base_size = base_size, base_family = "Arial") +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "#D9DEE2", linewidth = 0.28),
      panel.border = element_rect(color = "#4D555A", linewidth = 0.45),
      strip.background = element_rect(
        fill = "#EEF1F3",
        color = "#9FAAB1",
        linewidth = 0.4
      ),
      strip.text.y.left = element_text(
        face = "bold",
        color = "#202326",
        angle = 0,
        size = 7.4
      ),
      axis.text = element_text(color = "#202326"),
      axis.text.y = element_text(size = 7.7),
      axis.title = element_text(color = "#202326", size = 8.5),
      plot.title = element_text(face = "bold", size = 10.5, color = "#202326"),
      plot.subtitle = element_text(size = 7.8, color = "#4F5960"),
      plot.title.position = "plot",
      legend.position = "bottom",
      legend.box = "vertical",
      legend.direction = "horizontal",
      legend.title = element_text(size = 7.5, face = "bold"),
      legend.text = element_text(size = 7),
      legend.key.height = grid::unit(3.5, "mm"),
      legend.key.width = grid::unit(5, "mm"),
      legend.spacing.x = grid::unit(2, "mm"),
      legend.margin = margin(t = -3, unit = "pt"),
      plot.margin = margin(7, 8, 6, 7, unit = "pt")
    )
}

plot_scales <- function(
    title,
    subtitle,
    scale_data,
    include_treatment_legend = FALSE) {
  significant_fdr <- scale_data$minus_log10_FDR[scale_data$FDR <= 0.25]
  if (!length(significant_fdr)) {
    stop("No terms at FDR <= 0.25 are available for the red color scale.")
  }
  fdr_limits <- range(significant_fdr, na.rm = TRUE)
  if (diff(fdr_limits) == 0) {
    fdr_limits <- fdr_limits + c(-0.05, 0.05)
  }
  fdr_breaks <- seq(fdr_limits[1], fdr_limits[2], length.out = 3)
  has_nonsignificant <- any(scale_data$FDR > 0.25)

  guides_list <- list(
    fill = guide_colorbar(
      order = 1,
      title.position = "top",
      barwidth = grid::unit(25, "mm"),
      barheight = grid::unit(2.8, "mm")
    ),
    size = guide_legend(
      order = 2,
      title.position = "top",
      nrow = 1,
      override.aes = list(shape = 21, fill = "#F4B3B0", color = "#202326")
    )
  )
  if (has_nonsignificant) {
    guides_list$colour <- guide_legend(
      order = 3,
      title.position = "top",
      nrow = 1,
      override.aes = list(
        shape = 21,
        fill = "#B8BDC3",
        colour = "#202326",
        size = 3.5
      )
    )
  }
  if (include_treatment_legend) {
    guides_list$shape <- guide_legend(
      order = 4,
      title.position = "top",
      nrow = 1,
      override.aes = list(size = 3.5, fill = "#C83E3A")
    )
  }

  list(
    geom_vline(
      xintercept = 0,
      color = "#4D555A",
      linewidth = 0.35,
      linetype = "dashed"
    ),
    scale_x_continuous(
      name = "Normalized enrichment score (NES)",
      limits = x_limits,
      breaks = x_breaks,
      expand = expansion(mult = c(0.03, 0.03))
    ),
    scale_size_continuous(
      name = "GO gene-set size",
      range = c(2.6, 6),
      limits = size_limits,
      breaks = size_breaks
    ),
    scale_fill_gradient(
      name = expression(-log[10] * "(FDR)"),
      low = "#FADBD8",
      high = "#B2182B",
      limits = fdr_limits,
      breaks = fdr_breaks,
      labels = scales::label_number(accuracy = 0.1)
    ),
    if (has_nonsignificant) {
      scale_color_manual(
        name = "FDR status",
        values = c("FDR > 0.25" = "#202326"),
        drop = TRUE
      )
    },
    facet_grid(
      rows = vars(Category),
      scales = "free_y",
      space = "free_y",
      switch = "y",
      drop = TRUE,
      labeller = as_labeller(
        c(
          "OXPHOS" = "OXPHOS",
          "Glycolysis / carbon catabolism" = "Glycolysis /\ncarbon catabolism",
          "ECM / fibrosis" = "ECM /\nfibrosis"
        )
      )
    ),
    labs(title = title, subtitle = subtitle, y = NULL),
    do.call(guides, guides_list),
    coord_cartesian(clip = "off"),
    common_theme()
  )
}

single_treatment_plot <- function(data, title) {
  significant <- data[data$FDR <= 0.25, , drop = FALSE]
  nonsignificant <- data[data$FDR > 0.25, , drop = FALSE]

  plot <- ggplot(
    data,
    aes(x = NES, y = Term_label, size = Set_size)
  ) +
    geom_point(
      data = significant,
      aes(fill = minus_log10_FDR),
      shape = 21,
      color = "#202326",
      stroke = 0.45,
      alpha = 0.96
    )

  if (nrow(nonsignificant)) {
    plot <- plot +
      geom_point(
        data = nonsignificant,
        aes(color = FDR_status),
        shape = 21,
        fill = "#B8BDC3",
        stroke = 0.45,
        alpha = 0.96
      )
  }

  plot +
    plot_scales(
      title = title,
      subtitle = paste0(
        "edgeR exact test; ComBat-seq PCS-adjusted counts\n",
        "red scale uses this plot's FDR range; gray = FDR > 0.25"
      ),
      scale_data = data
    )
}

figure_2dg <- single_treatment_plot(plot_2dg, "2DG + OX vs control")
figure_ldha <- single_treatment_plot(plot_ldha, "LDHA + OX vs control")

figure_combined <- ggplot(
  plot_combined,
  aes(
    x = NES,
    y = Term_label,
    size = Set_size,
    shape = Treatment,
    group = Treatment
  )
) +
  geom_point(
    data = plot_combined[
      plot_combined$Treatment == "2DG + OX" &
        plot_combined$FDR <= 0.25,
      ,
      drop = FALSE
    ],
    aes(fill = minus_log10_FDR),
    position = position_nudge(y = 0.13),
    color = "#202326",
    stroke = 0.45,
    alpha = 0.96
  ) +
  geom_point(
    data = plot_combined[
      plot_combined$Treatment == "LDHA + OX" &
        plot_combined$FDR <= 0.25,
      ,
      drop = FALSE
    ],
    aes(fill = minus_log10_FDR),
    position = position_nudge(y = -0.13),
    color = "#202326",
    stroke = 0.45,
    alpha = 0.96
  ) +
  geom_point(
    data = plot_combined[
      plot_combined$Treatment == "2DG + OX" &
        plot_combined$FDR > 0.25,
      ,
      drop = FALSE
    ],
    aes(color = FDR_status),
    position = position_nudge(y = 0.13),
    fill = "#B8BDC3",
    stroke = 0.45,
    alpha = 0.96
  ) +
  geom_point(
    data = plot_combined[
      plot_combined$Treatment == "LDHA + OX" &
        plot_combined$FDR > 0.25,
      ,
      drop = FALSE
    ],
    aes(color = FDR_status),
    position = position_nudge(y = -0.13),
    fill = "#B8BDC3",
    stroke = 0.45,
    alpha = 0.96
  ) +
  scale_shape_manual(
    name = "Treatment",
    values = c("2DG + OX" = 21, "LDHA + OX" = 24)
  ) +
  plot_scales(
    title = "2DG + OX and LDHA + OX vs control",
    subtitle = paste0(
      "edgeR exact test; union of both term lists\nred scale uses this ",
      "plot's FDR range; gray = FDR > 0.25"
    ),
    scale_data = plot_combined,
    include_treatment_legend = TRUE
  )

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
    figure_2dg,
    "2DG_OX_edgeR_exact_selected_GSEA_dotplot",
    width_mm = 128,
    height_mm = 155
  ),
  save_plot(
    figure_ldha,
    "LDHA_OX_edgeR_exact_selected_GSEA_dotplot",
    width_mm = 128,
    height_mm = 150
  ),
  save_plot(
    figure_combined,
    "2DG_OX_and_LDHA_OX_edgeR_exact_union_GSEA_dotplot",
    width_mm = 152,
    height_mm = 170
  )
)

export_columns <- c(
  "Treatment", "Category", "GO_ID", "Term", "Ontology", "Direction",
  "NES", "Enrichment_score", "Set_size", "Core_gene_count",
  "P_value", "FDR", "q_value", "FDR_status", "Core_genes"
)
export_2dg <- plot_2dg[order(plot_2dg$Treatment_order), export_columns]
export_ldha <- plot_ldha[order(plot_ldha$Treatment_order), export_columns]
export_combined <- plot_combined[
  order(plot_combined$Category, plot_combined$Term_label, plot_combined$Treatment),
  export_columns
]

csv_2dg <- file.path(output_dir, "2DG_OX_edgeR_exact_selected_GSEA_values.csv")
csv_ldha <- file.path(output_dir, "LDHA_OX_edgeR_exact_selected_GSEA_values.csv")
csv_combined <- file.path(
  output_dir,
  "2DG_OX_and_LDHA_OX_edgeR_exact_union_GSEA_values.csv"
)
xlsx_file <- file.path(
  output_dir,
  "2DG_OX_and_LDHA_OX_edgeR_exact_selected_GSEA_values.xlsx"
)

write.csv(export_2dg, csv_2dg, row.names = FALSE, na = "")
write.csv(export_ldha, csv_ldha, row.names = FALSE, na = "")
write.csv(export_combined, csv_combined, row.names = FALSE, na = "")
write.xlsx(
  list(
    `2DG_edgeR_exact` = export_2dg,
    `LDHA_edgeR_exact` = export_ldha,
    `combined_union` = export_combined
  ),
  xlsx_file,
  asTable = TRUE,
  overwrite = TRUE
)

manifest_file <- file.path(output_dir, "figure_manifest.txt")
writeLines(
  c(
    "Figures: selected GO GSEA terms for 2DG + OX and LDHA + OX",
    "Differential-expression method: edgeR exact test only",
    "Counts: ComBat-seq PCS-adjusted",
    "Ranking: sign(log2FC) x -log10(raw P value)",
    "X axis: normalized enrichment score (NES)",
    "Dot size: GO gene-set size",
    "Red dot fill: -log10(FDR) for FDR <= 0.25",
    "Gray dot fill: FDR > 0.25 only",
    "Color scale: independent significant-term min/max for each figure",
    "Glycolysis and carbon catabolism are combined into one group",
    "Removed term: glycolytic process through glucose-6-phosphate",
    "Combined plot: full union of terms from both treatment-specific lists",
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    "2DG + OX terms:",
    paste0("- ", groups_2DG_IACS, ": ", terms_2DG_IACS),
    "",
    "LDHA + OX terms:",
    paste0("- ", groups_GNE_IACS, ": ", terms_GNE_IACS)
  ),
  manifest_file
)

cat(
  "Created:\n",
  paste(
    c(
      outputs,
      csv_2dg,
      csv_ldha,
      csv_combined,
      xlsx_file,
      manifest_file
    ),
    collapse = "\n"
  ),
  "\n"
)
