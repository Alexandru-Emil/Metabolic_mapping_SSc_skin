suppressPackageStartupMessages({
  library(SpatialFeatureExperiment)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(patchwork)
  library(ggsankey)
  library(pheatmap)
  library(viridisLite)
  library(scales)
  library(grid)
  library(svglite)
})

set.seed(20260729)
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 1:2) {
  stop(
    "Usage: 05_generate_seven_donor_manuscript_panels.R ",
    "<spe2025-12-17 umap.RData> [report_dir]"
  )
}

script_argument <- grep("^--file=", commandArgs(), value = TRUE)
stopifnot(length(script_argument) == 1L)
project_root <- normalizePath(
  if (length(args) == 2L) {
    args[[2]]
  } else {
    dirname(sub("^--file=", "", script_argument))
  },
  winslash = "/",
  mustWork = TRUE
)
output_root <- file.path(project_root, "outputs")
score_dir <- file.path(output_root, "seven_donor_scores")
panel_dir <- file.path(output_root, "seven_donor_manuscript_panels")
table_dir <- file.path(panel_dir, "tables")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

sfe_path <- file.path(output_root, "sfe_seven_donors_annotated_final.rds")
energy_path <- file.path(
  score_dir,
  "seven_donor_scmetabolism_energy_cell_scores.csv"
)
ecm_path <- file.path(
  score_dir,
  "seven_donor_matrisome_aucell_cell_scores.csv"
)
major_umap_path <- file.path(output_root, "major_joint_labels_and_umap.csv")
fib_umap_path <- file.path(output_root, "fibroblast_joint_labels_and_umap.csv")
ec_umap_path <- file.path(output_root, "endothelial_joint_labels_and_umap.csv")
spe_path <- normalizePath(
  args[[1]],
  winslash = "/",
  mustWork = TRUE
)

sample_to_donor <- c(
  "Validation1" = "Validation1",
  "Validation2" = "Validation2",
  "ExcludedValidation" = "ExcludedValidation",
  "Validation3" = "Validation3",
  "Validation4" = "Validation4",
  "Validation5" = "Validation5",
  "Validation6" = "Validation6"
)
donor_levels <- unname(sample_to_donor)

major_levels <- c(
  "epithelial",
  "fibroblast",
  "pericyte",
  "endothelial",
  "muscle",
  "telocyte",
  "Schwann cell",
  "myeloid cell",
  "lymphocyte"
)
fib_levels <- c(
  "papillary_Fib",
  "COMP_Fib",
  "COL8A1_Fib",
  "PI16_Fib",
  "CXCL12_Fib",
  "CCL19_Fib",
  "NGFR_Fib",
  "COCH_Fib"
)
ec_levels <- c(
  "cycling_EC",
  "ACKR1_EC",
  "EPC",
  "ACTA2_EC",
  "HEY1_EC",
  "lymphatic_EC"
)
fib_met_levels <- c("Met_hi_Fib", "Other_Fib")
ec_met_levels <- c("Met_hi_EC", "Other_EC")

fib_colors <- c("Met_hi_Fib" = "#F8766D", "Other_Fib" = "#00BFC4")
ec_met_colors <- c("Met_hi_EC" = "#E69F00", "Other_EC" = "#56B4E9")
major_colors <- setNames(scales::hue_pal()(length(major_levels)), major_levels)
fib_subtype_colors <- setNames(scales::hue_pal()(length(fib_levels)), fib_levels)
ec_subtype_colors <- c(
  "cycling_EC" = "#E69F00",
  "ACKR1_EC" = "#D55E00",
  "EPC" = "#009E73",
  "ACTA2_EC" = "#CC79A7",
  "HEY1_EC" = "#0072B2",
  "lymphatic_EC" = "#999999"
)
score_palette <- c("#2C7BB6", "#ABD9E9", "#FFFFBF", "#FDAE61", "#D7191C")
heatmap_palette <- viridisLite::viridis(100)

glycolysis_markers <- c("GLUT1", "HK1", "PFKL1_PFKM", "PKM2", "LDHA", "GAPDH")
tca_markers <- c("CS", "OGDH", "ATP5A", "SDHA")
protein_markers <- c("aSMA", "Cdh11", "FAP")
metabolic_markers <- c(
  glycolysis_markers,
  tca_markers,
  "CD98", "ACAC", "CPT1a", "G6PD", "pAMPK", "p_mTOR",
  "NOX4", "HIF1a", "PGC1a", "mtTFA", "pNRF2"
)
major_markers <- c(
  "KRT5", "CDH1", "COL1A1", "PDGFRA", "RGS5", "PDGFRB",
  "VWF", "ENG", "MYL9", "DES", "FOXL1", "PLN", "MPZ", "SCN7A",
  "MRC1", "CD1C", "CD2", "MS4A1"
)
fib_gene_markers <- c(
  "ACTA2", "COL8A1", "FN1", "SFRP4", "CCL19", "APOE",
  "COL18A1", "WIF1", "IL6", "IL10", "CCL20", "IL1R1",
  "TNC", "LOX", "LOXL2", "MMP9", "LUM", "PDGFRA",
  "PRG4", "PI16", "CD34", "FBLN2"
)
ec_gene_markers <- c(
  "MKI67", "TOP2A", "ACKR1", "SELP", "KDR", "APLNR",
  "FAP", "ACTA2", "HEY1", "GJA5", "LYVE1", "PDPN"
)

fib_imc_umap_name <- "umap_n_neighbors_15_cell_type_Stromal_1"
ec_imc_umap_name <- "umap_n_neighbors_15_cell_type_Endothelial_3"

panel_manifest <- tribble(
  ~panel, ~file_stem, ~width_in, ~height_in, ~color_contract,
  "FibMET UMAP", "fibmet_umap_seven_donors", 4.2, 3.6, "salmon/cyan",
  "FibMET IMC score UMAPs", "fibmet_imc_score_umaps_seven_donors", 6.8, 4.0, "blue-yellow-red [-1,1]",
  "FibMET protein expression", "fibmet_protein_expression_seven_donors", 6.8, 3.2, "salmon/cyan",
  "FibMET IMC heatmap", "fibmet_imc_heatmap_seven_donors", 4.8, 3.35, "viridis [-1,0.5]",
  "FibMET sample scores", "fibmet_sample_scores_seven_donors", 5.2, 3.7, "salmon/cyan",
  "FibMET Xenium gene dot plot", "fibmet_xenium_gene_dotplot_seven_donors", 4.8, 4.4, "blue-white-red [-1,1]",
  "FibMET Sankey", "fibmet_xenium_to_imc_sankey_seven_donors", 5.0, 3.7, "subtype categorical",
  "FibMET ECM dot plot", "fibmet_ecm_dotplot_seven_donors", 4.2, 4.2, "white-red",
  "FibMET scMetabolism", "fibmet_scmetabolism_dotplot_seven_donors", 5.6, 3.4, "white-red",
  "FibMET IMC-Xenium correlation", "fibmet_imc_scmetabolism_correlation_seven_donors", 7.2, 4.8, "salmon/cyan",
  "EndMET UMAP", "endmet_umap_seven_donors", 4.2, 3.6, "orange/blue",
  "EndMET IMC heatmap", "endmet_imc_heatmap_seven_donors", 4.8, 3.35, "viridis [-1,1]",
  "EndMET scMetabolism", "endmet_scmetabolism_dotplot_seven_donors", 5.6, 3.4, "white-red",
  "Xenium EC UMAP", "xenium_ec_umap_seven_donors", 4.2, 3.6, "EC categorical",
  "Xenium EC gene dot plot", "xenium_ec_gene_dotplot_seven_donors", 5.4, 3.8, "white-red",
  "Xenium EC ECM violins", "xenium_ec_ecm_violin_seven_donors", 6.5, 4.6, "EC categorical",
  "ACTA2_EC-EndMET Sankey", "acta2_ec_to_endmet_sankey_seven_donors", 4.8, 3.4, "pink/orange/blue",
  "ACTA2_EC scMetabolism", "acta2_ec_scmetabolism_dotplot_seven_donors", 5.6, 3.4, "white-red",
  "EndMET IMC-Xenium correlation", "endmet_imc_scmetabolism_correlation_seven_donors", 7.2, 4.8, "orange/blue",
  "Xenium major UMAP", "xenium_major_umap_seven_donors", 6.2, 4.4, "major categorical",
  "Xenium major marker dot plot", "xenium_major_marker_dotplot_seven_donors", 5.8, 4.4, "white-red"
)
write.csv(
  panel_manifest,
  file.path(panel_dir, "seven_donor_figure_export_manifest.csv"),
  row.names = FALSE
)

save_plot_set <- function(plot, stem, width, height) {
  ggsave(
    file.path(panel_dir, paste0(stem, ".png")),
    plot,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )
  ggsave(
    file.path(panel_dir, paste0(stem, ".pdf")),
    plot,
    width = width,
    height = height,
    device = cairo_pdf,
    bg = "white"
  )
  ggsave(
    file.path(panel_dir, paste0(stem, ".svg")),
    plot,
    width = width,
    height = height,
    device = svglite::svglite,
    bg = "white"
  )
}

save_grob_set <- function(grob, stem, width, height) {
  devices <- list(
    png = function(path) png(path, width, height, units = "in", res = 300, bg = "white"),
    pdf = function(path) cairo_pdf(path, width, height, bg = "white"),
    svg = function(path) svglite::svglite(path, width, height, bg = "white")
  )
  for (extension in names(devices)) {
    path <- file.path(panel_dir, paste0(stem, ".", extension))
    devices[[extension]](path)
    grid.newpage()
    grid.draw(grob)
    dev.off()
  }
}

rescale01 <- function(x) {
  limits <- range(x, na.rm = TRUE)
  if (!all(is.finite(limits)) || diff(limits) == 0) return(rep(0.5, length(x)))
  (x - limits[1]) / diff(limits)
}

format_p <- function(p) {
  ifelse(
    is.na(p),
    "p = NA",
    ifelse(p < 0.001, "p < 0.001", paste0("p = ", formatC(p, digits = 3, format = "f")))
  )
}

base_umap_theme <- theme_classic(base_size = 8) +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6.8, color = "grey20"),
    axis.ticks = element_line(linewidth = 0.35),
    axis.line = element_line(linewidth = 0.4),
    legend.title = element_text(size = 6.6),
    legend.text = element_text(size = 6.8),
    legend.key.height = unit(0.25, "cm"),
    legend.key.width = unit(0.28, "cm"),
    legend.spacing.y = unit(0.02, "cm"),
    legend.margin = margin(0, 0, 0, 1),
    plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
    plot.margin = margin(3, 4, 3, 3)
  )

paired_wilcox <- function(data, feature, value, group, group1, group2) {
  bind_rows(lapply(unique(data[[feature]]), function(current_feature) {
    wide <- data %>%
      filter(.data[[feature]] == current_feature) %>%
      select(donor, all_of(group), all_of(value)) %>%
      pivot_wider(names_from = all_of(group), values_from = all_of(value)) %>%
      filter(is.finite(.data[[group1]]), is.finite(.data[[group2]]))
    test <- if (nrow(wide) >= 2L) {
      suppressWarnings(wilcox.test(wide[[group1]], wide[[group2]], paired = TRUE, exact = FALSE))
    } else {
      NULL
    }
    tibble(
      feature = current_feature,
      n_paired_donors = nrow(wide),
      p_value = if (is.null(test)) NA_real_ else test$p.value
    )
  }))
}

dot_summary <- function(sfe, cells, group_values, markers) {
  markers <- markers[markers %in% rownames(sfe)]
  groups <- unique(as.character(group_values))
  groups <- groups[!is.na(groups)]
  expression <- assay(sfe, "logcounts")[markers, cells, drop = FALSE]
  counts <- assay(sfe, "counts")[markers, cells, drop = FALSE]
  bind_rows(lapply(groups, function(current_group) {
    idx <- as.character(group_values) == current_group
    tibble(
      group = current_group,
      gene = markers,
      mean_expression = Matrix::rowMeans(expression[, idx, drop = FALSE]),
      percent_expressed = 100 * Matrix::rowMeans(counts[, idx, drop = FALSE] > 0),
      n_cells = sum(idx)
    )
  })) %>%
    group_by(gene) %>%
    mutate(
      scaled_expression = if (n() > 1L && sd(mean_expression) > 0) {
        as.numeric(base::scale(mean_expression))
      } else {
        0
      }
    ) %>%
    ungroup()
}

make_energy_dotplot <- function(summary, group_col, group_levels, title) {
  plot_data <- summary %>%
    mutate(
      group = factor(.data[[group_col]], levels = group_levels),
      pathway = factor(pathway, levels = rev(c("Glycolysis", "TCA / OXPHOS"))),
      size_score = rescale01(mean_score)
    )
  ggplot(plot_data, aes(x = group, y = pathway)) +
    geom_point(
      aes(fill = mean_score, size = size_score),
      shape = 21,
      color = "grey35",
      stroke = 0.35
    ) +
    scale_fill_gradient(low = "white", high = "#D7301F", name = "Mean AUCell score") +
    scale_size(range = c(6.5, 12.5), guide = "none") +
    labs(title = title, x = NULL, y = NULL) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(size = 13, face = "bold"),
      axis.text.y = element_text(size = 13),
      legend.title = element_text(size = 11.5),
      legend.text = element_text(size = 10.5),
      panel.grid.major = element_line(color = "grey88"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 14.5, hjust = 0.5),
      plot.margin = margin(5.5, 6, 5.5, 5.5)
    )
}

make_sankey <- function(data, source_col, target_col, palette, title) {
  sankey_data <- data %>%
    transmute(
      source = as.character(.data[[source_col]]),
      target = as.character(.data[[target_col]])
    ) %>%
    filter(!is.na(source), !is.na(target)) %>%
    ggsankey::make_long(source, target)

  ggplot(
    sankey_data,
    aes(
      x = x,
      next_x = next_x,
      node = node,
      next_node = next_node,
      fill = factor(node),
      label = node
    )
  ) +
    geom_sankey(flow.alpha = 0.42, node.color = "grey35", width = 0.16) +
    geom_sankey_label(size = 2.7, color = "grey15", fill = "white") +
    scale_fill_manual(values = palette, guide = "none") +
    scale_x_discrete(labels = c("Xenium", "IMC"), expand = c(0.08, 0.08)) +
    labs(title = title, x = NULL, y = NULL) +
    theme_sankey(base_size = 9) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 8, face = "bold"),
      plot.margin = margin(4, 4, 4, 4)
    )
}

make_correlation <- function(data, group_col, group_levels, colors, title) {
  stats <- data %>%
    group_by(pathway, .data[[group_col]]) %>%
    summarise(
      n = n(),
      rho = if (n() >= 3L) {
        suppressWarnings(cor(imc_score, xenium_score, method = "spearman"))
      } else {
        NA_real_
      },
      .groups = "drop"
    ) %>%
    mutate(label = paste0("n=", n, ", rho=", sprintf("%.2f", rho)))

  label_positions <- data %>%
    group_by(pathway, .data[[group_col]]) %>%
    summarise(
      x = quantile(imc_score, 0.04, na.rm = TRUE),
      y = quantile(xenium_score, 0.96, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(stats, by = c("pathway", group_col))

  plot <- ggplot(
    data,
    aes(
      x = imc_score,
      y = xenium_score,
      color = .data[[group_col]]
    )
  ) +
    geom_point(size = 0.55, alpha = 0.32, stroke = 0) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.65) +
    geom_text(
      data = label_positions,
      aes(x = x, y = y, label = label, color = .data[[group_col]]),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 2.7
    ) +
    facet_grid(rows = vars(.data[[group_col]]), cols = vars(pathway), scales = "free") +
    scale_color_manual(values = colors, limits = group_levels, drop = FALSE) +
    labs(
      title = title,
      subtitle = "One point per cell with a validated shared IMC-Xenium identifier",
      x = "IMC protein score",
      y = "Xenium scMetabolism AUCell score",
      color = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 9, face = "bold"),
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 8.5, hjust = 0.5),
      panel.grid.minor = element_blank()
    )
  list(plot = plot, stats = stats)
}

message("Reading seven-donor SFE, score caches, and IMC SPE")
sfe7 <- readRDS(sfe_path)
sfe_meta <- as.data.frame(colData(sfe7), check.names = FALSE)
sfe_meta$cell <- colnames(sfe7)

energy <- read.csv(energy_path, check.names = FALSE)
ecm <- read.csv(ecm_path, check.names = FALSE)
major_umap <- read.csv(major_umap_path, check.names = FALSE)
fib_joint_umap <- read.csv(fib_umap_path, check.names = FALSE)
ec_joint_umap <- read.csv(ec_umap_path, check.names = FALSE)

imc_env <- new.env(parent = emptyenv())
load(spe_path, envir = imc_env)
spe <- imc_env$spe
keep_imc <- as.character(spe$sample_id) %in% names(sample_to_donor)
spe7 <- spe[, keep_imc]
spe7$donor <- factor(
  unname(sample_to_donor[as.character(spe7$sample_id)]),
  levels = donor_levels
)
spe7$metfiblabel <- case_when(
  as.character(spe7$Rphenograph_500) %in% c("1", "2") ~ "Met_hi_Fib",
  as.character(spe7$Rphenograph_500) == "3" ~ "Other_Fib",
  TRUE ~ NA_character_
)
spe7$metEClabel <- case_when(
  as.character(spe7$Rphenograph_200EC) %in% c("1", "2") ~ "Met_hi_EC",
  as.character(spe7$Rphenograph_200EC) == "3" ~ "Other_EC",
  TRUE ~ NA_character_
)

stopifnot(
  length(unique(as.character(spe7$sample_id))) == 7L,
  all(c(fib_imc_umap_name, ec_imc_umap_name) %in% reducedDimNames(spe7)),
  all(unique(c(metabolic_markers, protein_markers)) %in% rownames(spe7))
)

asinh7 <- as.matrix(assay(spe7, "asinh"))
row_mean <- rowMeans(asinh7, na.rm = TRUE)
row_sd <- apply(asinh7, 1, sd, na.rm = TRUE)
row_sd[!is.finite(row_sd) | row_sd == 0] <- 1
z7 <- sweep(sweep(asinh7, 1, row_mean, "-"), 1, row_sd, "/")

glycolysis_score <- colMeans(z7[glycolysis_markers, , drop = FALSE])
tca_score <- colMeans(z7[tca_markers, , drop = FALSE])
names(glycolysis_score) <- colnames(spe7)
names(tca_score) <- colnames(spe7)

message("Generating Xenium-only annotation panels")
major_plot_data <- major_umap %>%
  select(cell, UMAP1, UMAP2) %>%
  left_join(
    sfe_meta %>% transmute(
      cell,
      donor = as.character(SampleId),
      annotation = as.character(lv1_anno_final)
    ),
    by = "cell"
  ) %>%
  filter(!is.na(annotation)) %>%
  mutate(annotation = factor(annotation, levels = major_levels))

p_major_umap <- ggplot(
  major_plot_data,
  aes(UMAP1, UMAP2, color = annotation)
) +
  geom_point(size = 0.22, alpha = 0.9, stroke = 0) +
  scale_color_manual(values = major_colors, drop = FALSE) +
  labs(x = "UMAP_1", y = "UMAP_2", color = NULL) +
  base_umap_theme +
  theme(
    legend.position = "right",
    legend.key.height = unit(0.28, "cm"),
    legend.text = element_text(size = 6.7)
  ) +
  guides(color = guide_legend(override.aes = list(size = 2.1, alpha = 1)))
save_plot_set(p_major_umap, "xenium_major_umap_seven_donors", 6.2, 4.4)
write.csv(
  major_plot_data,
  file.path(table_dir, "xenium_major_umap_seven_donors.csv"),
  row.names = FALSE
)

major_cells <- which(!is.na(sfe7$lv1_anno_final))
major_dot <- dot_summary(
  sfe7,
  major_cells,
  sfe7$lv1_anno_final[major_cells],
  major_markers
) %>%
  mutate(
    group = factor(group, levels = rev(major_levels)),
    gene = factor(gene, levels = major_markers)
  )
major_color_limit <- max(1.6, quantile(major_dot$mean_expression, 0.98, na.rm = TRUE))
p_major_dot <- ggplot(
  major_dot,
  aes(gene, group, size = percent_expressed, color = mean_expression)
) +
  geom_point() +
  scale_size(
    name = "Percent Expressed",
    range = c(1.2, 5.2),
    limits = c(0, 100),
    breaks = c(25, 50, 75, 100)
  ) +
  scale_color_gradient(
    low = "#F7F7F7",
    high = "#A50000",
    limits = c(0, major_color_limit),
    oob = squish,
    name = "Average\nExpression"
  ) +
  labs(x = "Features", y = "Identity") +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_text(size = 6.4, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6.9),
    axis.title = element_text(size = 7.5),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6.5),
    plot.margin = margin(3, 4, 4, 3)
  )
save_plot_set(p_major_dot, "xenium_major_marker_dotplot_seven_donors", 5.8, 4.4)
write.csv(
  major_dot,
  file.path(table_dir, "xenium_major_marker_dotplot_seven_donors.csv"),
  row.names = FALSE
)

message("Generating seven-donor FibMET IMC panels")
fib_embedding <- reducedDim(spe7, fib_imc_umap_name)
fib_imc <- tibble(
  cell = colnames(spe7),
  UMAP1 = fib_embedding[, 1],
  UMAP2 = fib_embedding[, 2],
  sample_id = as.character(spe7$sample_id),
  donor = as.character(spe7$donor),
  metfiblabel = as.character(spe7$metfiblabel),
  glycolysis_score = unname(glycolysis_score[colnames(spe7)]),
  tca_score = unname(tca_score[colnames(spe7)])
) %>%
  filter(
    complete.cases(UMAP1, UMAP2),
    metfiblabel %in% fib_met_levels
  ) %>%
  mutate(
    donor = factor(donor, levels = donor_levels),
    metfiblabel = factor(metfiblabel, levels = fib_met_levels)
  )

p_fib_umap <- ggplot(
  fib_imc,
  aes(UMAP1, UMAP2, color = metfiblabel)
) +
  geom_point(size = 0.55, alpha = 0.9, stroke = 0) +
  scale_color_manual(values = fib_colors, drop = FALSE) +
  labs(
    title = "Metabolic phenotypes of fibroblasts (FibMET)",
    x = "UMAP1",
    y = "UMAP2",
    color = NULL
  ) +
  base_umap_theme +
  guides(color = guide_legend(override.aes = list(size = 2.2, alpha = 1)))
save_plot_set(p_fib_umap, "fibmet_umap_seven_donors", 4.2, 3.6)

make_score_umap <- function(data, score, title) {
  ggplot(data, aes(UMAP1, UMAP2, color = .data[[score]])) +
    geom_point(size = 0.52, alpha = 0.92, stroke = 0) +
    scale_color_gradientn(
      colors = score_palette,
      limits = c(-1, 1),
      breaks = seq(-1, 1, 0.5),
      oob = squish,
      name = "Score",
      guide = guide_colorbar(
        barwidth = unit(0.18, "cm"),
        barheight = unit(0.95, "cm"),
        title.position = "top"
      )
    ) +
    labs(title = title, x = "UMAP1", y = "UMAP2") +
    base_umap_theme +
    theme(
      legend.title = element_text(size = 6.3),
      legend.text = element_text(size = 5.9),
      plot.margin = margin(3, 2, 3, 2)
    )
}
p_fib_scores <- (
  make_score_umap(fib_imc, "glycolysis_score", "Glycolysis score") +
    make_score_umap(fib_imc, "tca_score", "TCA/OXPHOS score")
) +
  plot_layout(guides = "keep") +
  plot_annotation(
    title = "Metabolic scores",
    theme = theme(
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5)
    )
  )
save_plot_set(p_fib_scores, "fibmet_imc_score_umaps_seven_donors", 6.8, 4.0)
write.csv(
  fib_imc,
  file.path(table_dir, "fibmet_imc_cells_seven_donors.csv"),
  row.names = FALSE
)

fib_idx <- match(fib_imc$cell, colnames(spe7))
fib_protein_cell <- as.data.frame(
  t(z7[protein_markers, fib_idx, drop = FALSE]),
  check.names = FALSE
) %>%
  rownames_to_column("cell") %>%
  left_join(fib_imc %>% select(cell, donor, metfiblabel), by = "cell") %>%
  pivot_longer(all_of(protein_markers), names_to = "marker", values_to = "expression")
fib_protein_sample <- fib_protein_cell %>%
  group_by(donor, metfiblabel, marker) %>%
  summarise(mean_expression = mean(expression), .groups = "drop") %>%
  mutate(
    donor = factor(donor, levels = donor_levels),
    metfiblabel = factor(metfiblabel, levels = fib_met_levels),
    marker = factor(marker, levels = protein_markers)
  )
fib_protein_stats <- paired_wilcox(
  fib_protein_sample,
  "marker",
  "mean_expression",
  "metfiblabel",
  "Met_hi_Fib",
  "Other_Fib"
)
fib_protein_pos <- fib_protein_sample %>%
  group_by(marker) %>%
  summarise(
    label_y = max(mean_expression) + 0.15 * max(diff(range(mean_expression)), 0.2),
    .groups = "drop"
  ) %>%
  left_join(
    fib_protein_stats %>% transmute(marker = feature, p_value),
    by = "marker"
  ) %>%
  mutate(label = format_p(p_value))
p_fib_protein <- ggplot(
  fib_protein_sample,
  aes(metfiblabel, mean_expression, fill = metfiblabel)
) +
  geom_violin(width = 0.78, alpha = 0.76, color = "grey25", linewidth = 0.28) +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", linewidth = 0.28) +
  geom_point(position = position_jitter(width = 0.045), size = 1.05, color = "black") +
  geom_text(
    data = fib_protein_pos,
    aes(x = 1.5, y = label_y, label = label),
    inherit.aes = FALSE,
    size = 2.05
  ) +
  facet_wrap(~marker, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = fib_colors) +
  labs(title = "Protein expression", x = NULL, y = "Mean expression z score") +
  theme_bw(base_size = 7.5) +
  theme(
    axis.text.x = element_text(size = 6.6, angle = 42, hjust = 1),
    axis.text.y = element_text(size = 6.2),
    strip.background = element_blank(),
    strip.text = element_text(size = 8, face = "bold"),
    legend.position = "none",
    plot.title = element_text(size = 9.5, face = "bold", hjust = 0.5),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.25)
  )
save_plot_set(p_fib_protein, "fibmet_protein_expression_seven_donors", 6.8, 3.2)
write.csv(
  fib_protein_sample,
  file.path(table_dir, "fibmet_protein_expression_sample_means_seven_donors.csv"),
  row.names = FALSE
)
write.csv(
  fib_protein_stats,
  file.path(table_dir, "fibmet_protein_expression_paired_wilcoxon_seven_donors.csv"),
  row.names = FALSE
)

make_imc_heatmap <- function(
    cells,
    group_values,
    group_levels,
    group_colors,
    breaks,
    stem) {
  cell_index <- match(cells, colnames(spe7))
  heat_data <- as.data.frame(
    t(z7[metabolic_markers, cell_index, drop = FALSE]),
    check.names = FALSE
  ) %>%
    mutate(group = as.character(group_values)) %>%
    pivot_longer(all_of(metabolic_markers), names_to = "marker", values_to = "value") %>%
    group_by(group, marker) %>%
    summarise(mean_expression = mean(value), .groups = "drop")
  matrix_df <- heat_data %>%
    pivot_wider(names_from = group, values_from = mean_expression)
  matrix_df <- matrix_df[match(metabolic_markers, matrix_df$marker), , drop = FALSE]
  matrix <- as.matrix(matrix_df[, group_levels, drop = FALSE])
  rownames(matrix) <- matrix_df$marker
  counts <- table(factor(group_values, levels = group_levels))
  annotation <- data.frame(
    ids = factor(group_levels, levels = group_levels),
    ncells = as.numeric(counts[group_levels]),
    row.names = group_levels
  )
  object <- pheatmap(
    matrix,
    color = heatmap_palette,
    breaks = seq(breaks[1], breaks[2], length.out = 101),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_col = annotation,
    annotation_colors = list(
      ids = group_colors,
      ncells = colorRampPalette(c("#F7F0F5", "#CC79A7"))(100)
    ),
    show_colnames = FALSE,
    border_color = NA,
    legend_breaks = pretty(breaks, n = 4),
    fontsize = 7,
    fontsize_row = 6.2,
    treeheight_row = 0,
    treeheight_col = 0,
    silent = TRUE
  )
  save_grob_set(object$gtable, stem, 4.8, 3.35)
  write.csv(
    heat_data,
    file.path(table_dir, paste0(stem, ".csv")),
    row.names = FALSE
  )
  object$gtable
}
fib_heatmap_grob <- make_imc_heatmap(
  fib_imc$cell,
  fib_imc$metfiblabel,
  fib_met_levels,
  fib_colors,
  c(-1, 0.5),
  "fibmet_imc_heatmap_seven_donors"
)

fib_score_sample <- fib_imc %>%
  select(donor, metfiblabel, Glycolysis = glycolysis_score, `TCA/OXPHOS` = tca_score) %>%
  pivot_longer(c(Glycolysis, `TCA/OXPHOS`), names_to = "pathway", values_to = "score") %>%
  group_by(donor, metfiblabel, pathway) %>%
  summarise(mean_score = mean(score), .groups = "drop") %>%
  mutate(
    metfiblabel = factor(metfiblabel, levels = fib_met_levels),
    pathway = factor(pathway, levels = c("Glycolysis", "TCA/OXPHOS"))
  )
fib_score_stats <- paired_wilcox(
  fib_score_sample,
  "pathway",
  "mean_score",
  "metfiblabel",
  "Met_hi_Fib",
  "Other_Fib"
)
make_score_box <- function(pathway_name) {
  data <- fib_score_sample %>% filter(pathway == pathway_name)
  stat <- fib_score_stats %>% filter(feature == pathway_name)
  y <- max(data$mean_score) + 0.15 * max(diff(range(data$mean_score)), 0.2)
  ggplot(data, aes(metfiblabel, mean_score, fill = metfiblabel)) +
    geom_boxplot(width = 0.52, outlier.shape = NA, alpha = 0.7, linewidth = 0.3) +
    geom_point(position = position_jitter(width = 0.045), size = 1.15, color = "black") +
    annotate("text", x = 1.5, y = y, label = format_p(stat$p_value), size = 2.15) +
    scale_fill_manual(values = fib_colors) +
    labs(title = pathway_name, x = NULL, y = "Mean IMC protein score") +
    theme_bw(base_size = 7.5) +
    theme(
      axis.text.x = element_text(size = 6.6, angle = 20, hjust = 1),
      legend.position = "none",
      plot.title = element_text(size = 8.5, face = "bold", hjust = 0.5),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.25)
    )
}
p_fib_sample <- make_score_box("Glycolysis") / make_score_box("TCA/OXPHOS")
save_plot_set(p_fib_sample, "fibmet_sample_scores_seven_donors", 5.2, 3.7)
write.csv(
  fib_score_sample,
  file.path(table_dir, "fibmet_sample_scores_seven_donors.csv"),
  row.names = FALSE
)
write.csv(
  fib_score_stats,
  file.path(table_dir, "fibmet_sample_score_tests_seven_donors.csv"),
  row.names = FALSE
)

message("Generating paired FibMET-Xenium panels")
fib_paired_cells <- which(
  as.character(sfe7$lv1_anno_final) == "fibroblast" &
    as.character(sfe7$metfiblabel) %in% fib_met_levels
)
fib_gene_dot <- dot_summary(
  sfe7,
  fib_paired_cells,
  sfe7$metfiblabel[fib_paired_cells],
  fib_gene_markers
) %>%
  mutate(
    group = factor(group, levels = rev(fib_met_levels)),
    gene = factor(gene, levels = fib_gene_markers)
  )
p_fib_gene_dot <- ggplot(
  fib_gene_dot,
  aes(group, gene, size = percent_expressed, color = scaled_expression)
) +
  geom_point() +
  scale_size(
    name = "Percent\nexpressed",
    range = c(1.4, 5),
    limits = c(0, 100),
    breaks = c(10, 20, 30, 40)
  ) +
  scale_color_gradient2(
    low = "#2C7BB6",
    mid = "white",
    high = "#D7191C",
    midpoint = 0,
    limits = c(-1, 1),
    oob = squish,
    name = "Average\nexpression"
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_text(size = 7, face = "bold"),
    axis.text.y = element_text(size = 6.7),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6.5)
  )
save_plot_set(p_fib_gene_dot, "fibmet_xenium_gene_dotplot_seven_donors", 4.8, 4.4)
write.csv(
  fib_gene_dot,
  file.path(table_dir, "fibmet_xenium_gene_dotplot_seven_donors.csv"),
  row.names = FALSE
)

fib_sankey_data <- sfe_meta %>%
  filter(
    lv1_anno_final == "fibroblast",
    lv2_anno_final %in% fib_levels,
    metfiblabel %in% fib_met_levels
  )
fib_sankey_palette <- c(
  fib_subtype_colors,
  "Met_hi_Fib" = "#F564C8",
  "Other_Fib" = "#F564C8"
)
p_fib_sankey <- make_sankey(
  fib_sankey_data,
  "lv2_anno_final",
  "metfiblabel",
  fib_sankey_palette,
  NULL
)
save_plot_set(p_fib_sankey, "fibmet_xenium_to_imc_sankey_seven_donors", 5.0, 3.7)
fib_sankey_counts <- fib_sankey_data %>%
  count(lv2_anno_final, metfiblabel, name = "Freq") %>%
  rename(Xenium = lv2_anno_final, IMC = metfiblabel)
write.csv(
  fib_sankey_counts,
  file.path(table_dir, "fibmet_xenium_to_imc_sankey_counts_seven_donors.csv"),
  row.names = FALSE
)

fib_ecm <- ecm %>%
  filter(
    ecm_score == "CoreMatrisome score",
    lv1_anno == "fibroblast",
    lv2_anno %in% fib_levels,
    metfiblabel %in% fib_met_levels
  ) %>%
  group_by(lv2_anno, metfiblabel) %>%
  summarise(mean_score = mean(score), n_cells = n(), .groups = "drop") %>%
  mutate(
    lv2_anno = factor(lv2_anno, levels = rev(fib_levels)),
    metfiblabel = factor(metfiblabel, levels = fib_met_levels),
    size_score = rescale01(mean_score)
  )
p_fib_ecm <- ggplot(
  fib_ecm,
  aes(metfiblabel, lv2_anno, size = size_score, color = mean_score)
) +
  geom_point() +
  scale_size(range = c(2, 7), guide = "none") +
  scale_color_gradient(low = "#F7F7F7", high = "#A50000", name = "average score") +
  labs(title = "ECM scores", x = NULL, y = NULL) +
  theme_classic(base_size = 8.5) +
  theme(
    axis.text.x = element_text(size = 7.5, face = "bold"),
    axis.text.y = element_text(size = 7),
    plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6.5)
  )
save_plot_set(p_fib_ecm, "fibmet_ecm_dotplot_seven_donors", 4.2, 4.2)
write.csv(
  fib_ecm,
  file.path(table_dir, "fibmet_ecm_dotplot_seven_donors.csv"),
  row.names = FALSE
)

fib_energy_long <- energy %>%
  filter(
    lv1_anno == "fibroblast",
    metfiblabel %in% fib_met_levels
  ) %>%
  select(
    cell,
    donor,
    metfiblabel,
    Glycolysis = glycolysis_scmetabolism,
    `TCA / OXPHOS` = tca_scmetabolism
  ) %>%
  pivot_longer(c(Glycolysis, `TCA / OXPHOS`), names_to = "pathway", values_to = "score")
fib_energy_summary <- fib_energy_long %>%
  group_by(pathway, metfiblabel) %>%
  summarise(mean_score = mean(score), n_cells = n(), .groups = "drop")
p_fib_scmet <- make_energy_dotplot(
  fib_energy_summary,
  "metfiblabel",
  c("Other_Fib", "Met_hi_Fib"),
  "Xenium-based scMetabolism scores"
)
save_plot_set(p_fib_scmet, "fibmet_scmetabolism_dotplot_seven_donors", 5.6, 3.4)
write.csv(
  fib_energy_summary,
  file.path(table_dir, "fibmet_scmetabolism_summary_seven_donors.csv"),
  row.names = FALSE
)

message("Generating seven-donor EndMET and Xenium EC panels")
ec_embedding <- reducedDim(spe7, ec_imc_umap_name)
ec_imc <- tibble(
  cell = colnames(spe7),
  UMAP1 = ec_embedding[, 1],
  UMAP2 = ec_embedding[, 2],
  sample_id = as.character(spe7$sample_id),
  donor = as.character(spe7$donor),
  metEClabel = as.character(spe7$metEClabel),
  glycolysis_score = unname(glycolysis_score[colnames(spe7)]),
  tca_score = unname(tca_score[colnames(spe7)])
) %>%
  filter(
    complete.cases(UMAP1, UMAP2),
    metEClabel %in% ec_met_levels
  ) %>%
  mutate(
    donor = factor(donor, levels = donor_levels),
    metEClabel = factor(metEClabel, levels = ec_met_levels)
  )
p_ec_imc_umap <- ggplot(ec_imc, aes(UMAP1, UMAP2, color = metEClabel)) +
  geom_point(size = 0.68, alpha = 0.92, stroke = 0) +
  scale_color_manual(values = ec_met_colors, drop = FALSE) +
  labs(
    title = "Metabolic phenotypes of\nendothelial cells (End-Met)",
    x = "UMAP_EC_1",
    y = "UMAP_EC_2",
    color = NULL
  ) +
  base_umap_theme +
  guides(color = guide_legend(override.aes = list(size = 2.2, alpha = 1)))
save_plot_set(p_ec_imc_umap, "endmet_umap_seven_donors", 4.2, 3.6)
write.csv(
  ec_imc,
  file.path(table_dir, "endmet_imc_cells_seven_donors.csv"),
  row.names = FALSE
)

ec_heatmap_grob <- make_imc_heatmap(
  ec_imc$cell,
  ec_imc$metEClabel,
  ec_met_levels,
  ec_met_colors,
  c(-1, 1),
  "endmet_imc_heatmap_seven_donors"
)

ec_energy_long <- energy %>%
  filter(
    lv1_anno == "endothelial",
    metEClabel %in% ec_met_levels
  ) %>%
  select(
    cell,
    donor,
    metEClabel,
    Glycolysis = glycolysis_scmetabolism,
    `TCA / OXPHOS` = tca_scmetabolism
  ) %>%
  pivot_longer(c(Glycolysis, `TCA / OXPHOS`), names_to = "pathway", values_to = "score")
ec_energy_summary <- ec_energy_long %>%
  group_by(pathway, metEClabel) %>%
  summarise(mean_score = mean(score), n_cells = n(), .groups = "drop")
p_ec_scmet <- make_energy_dotplot(
  ec_energy_summary,
  "metEClabel",
  c("Other_EC", "Met_hi_EC"),
  "Xenium-based scMetabolism scores"
)
save_plot_set(p_ec_scmet, "endmet_scmetabolism_dotplot_seven_donors", 5.6, 3.4)
write.csv(
  ec_energy_summary,
  file.path(table_dir, "endmet_scmetabolism_summary_seven_donors.csv"),
  row.names = FALSE
)

ec_plot_data <- ec_joint_umap %>%
  select(cell, UMAP1, UMAP2) %>%
  left_join(
    sfe_meta %>% transmute(
      cell,
      donor = as.character(SampleId),
      lv2_anno = as.character(lv2_anno_final)
    ),
    by = "cell"
  ) %>%
  filter(lv2_anno %in% ec_levels) %>%
  mutate(lv2_anno = factor(lv2_anno, levels = ec_levels))
p_ec_xenium_umap <- ggplot(ec_plot_data, aes(UMAP1, UMAP2, color = lv2_anno)) +
  geom_point(size = 0.34, alpha = 0.9, stroke = 0) +
  scale_color_manual(values = ec_subtype_colors, drop = FALSE) +
  labs(x = "UMAP_1", y = "UMAP_2", color = NULL) +
  base_umap_theme +
  theme(legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 2.1, alpha = 1)))
save_plot_set(p_ec_xenium_umap, "xenium_ec_umap_seven_donors", 4.2, 3.6)
write.csv(
  ec_plot_data,
  file.path(table_dir, "xenium_ec_umap_seven_donors.csv"),
  row.names = FALSE
)

ec_cells <- which(as.character(sfe7$lv2_anno_final) %in% ec_levels)
ec_gene_dot <- dot_summary(
  sfe7,
  ec_cells,
  sfe7$lv2_anno_final[ec_cells],
  ec_gene_markers
) %>%
  mutate(
    group = factor(group, levels = rev(ec_levels)),
    gene = factor(gene, levels = ec_gene_markers)
  )
ec_color_limit <- max(0.8, quantile(ec_gene_dot$mean_expression, 0.98, na.rm = TRUE))
p_ec_gene_dot <- ggplot(
  ec_gene_dot,
  aes(gene, group, size = percent_expressed, color = mean_expression)
) +
  geom_point() +
  scale_size(
    name = "Percent\nexpressed",
    range = c(1.2, 5),
    limits = c(0, 100),
    breaks = c(25, 50, 75, 100)
  ) +
  scale_color_gradient(
    low = "#F7F7F7",
    high = "#A50000",
    limits = c(0, ec_color_limit),
    oob = squish,
    name = "Average\nexpression"
  ) +
  labs(x = "Features", y = "Identity") +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_text(size = 6.4, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6.8),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6.5)
  )
save_plot_set(p_ec_gene_dot, "xenium_ec_gene_dotplot_seven_donors", 5.4, 3.8)
write.csv(
  ec_gene_dot,
  file.path(table_dir, "xenium_ec_gene_dotplot_seven_donors.csv"),
  row.names = FALSE
)

ec_ecm <- ecm %>%
  filter(
    lv1_anno == "endothelial",
    lv2_anno %in% ec_levels,
    ecm_score %in% c("CoreMatrisome score", "Collagens score")
  ) %>%
  mutate(
    lv2_anno = factor(lv2_anno, levels = ec_levels),
    ecm_score = factor(
      ecm_score,
      levels = c("CoreMatrisome score", "Collagens score"),
      labels = c("Core\nmatrisome", "Collagens")
    )
  )
p_ec_ecm <- ggplot(ec_ecm, aes(lv2_anno, score, fill = lv2_anno)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.72, color = "grey30", linewidth = 0.25) +
  geom_boxplot(width = 0.11, outlier.shape = NA, fill = "white", linewidth = 0.22) +
  facet_wrap(~ecm_score, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = ec_subtype_colors, drop = FALSE) +
  labs(x = NULL, y = "AUCell score", fill = NULL) +
  theme_bw(base_size = 8) +
  theme(
    axis.text.x = element_text(size = 6.4, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6.4),
    strip.background = element_blank(),
    strip.text = element_text(size = 7.4, face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )
save_plot_set(p_ec_ecm, "xenium_ec_ecm_violin_seven_donors", 6.5, 4.6)

acta2_sankey_data <- sfe_meta %>%
  filter(
    lv2_anno_final == "ACTA2_EC",
    metEClabel %in% ec_met_levels
  ) %>%
  mutate(acta2_group = "ACTA2_EC")
p_acta2_sankey <- make_sankey(
  acta2_sankey_data,
  "acta2_group",
  "metEClabel",
  c("ACTA2_EC" = "#E89A9A", ec_met_colors),
  "ACTA2_EC to End-Met state"
)
save_plot_set(p_acta2_sankey, "acta2_ec_to_endmet_sankey_seven_donors", 4.8, 3.4)
acta2_sankey_counts <- acta2_sankey_data %>%
  count(acta2_group, metEClabel, name = "Freq") %>%
  rename(Xenium = acta2_group, IMC = metEClabel)
write.csv(
  acta2_sankey_counts,
  file.path(table_dir, "acta2_ec_to_endmet_sankey_counts_seven_donors.csv"),
  row.names = FALSE
)

acta2_energy_long <- energy %>%
  filter(lv1_anno == "endothelial", lv2_anno %in% ec_levels) %>%
  mutate(acta2_group = ifelse(lv2_anno == "ACTA2_EC", "ACTA2_EC", "Other_EC")) %>%
  select(
    cell,
    donor,
    acta2_group,
    Glycolysis = glycolysis_scmetabolism,
    `TCA / OXPHOS` = tca_scmetabolism
  ) %>%
  pivot_longer(c(Glycolysis, `TCA / OXPHOS`), names_to = "pathway", values_to = "score")
acta2_energy_summary <- acta2_energy_long %>%
  group_by(pathway, acta2_group) %>%
  summarise(mean_score = mean(score), n_cells = n(), .groups = "drop")
p_acta2_scmet <- make_energy_dotplot(
  acta2_energy_summary,
  "acta2_group",
  c("Other_EC", "ACTA2_EC"),
  "ACTA2_EC scMetabolism scores"
)
save_plot_set(p_acta2_scmet, "acta2_ec_scmetabolism_dotplot_seven_donors", 5.6, 3.4)
write.csv(
  acta2_energy_summary,
  file.path(table_dir, "acta2_ec_scmetabolism_summary_seven_donors.csv"),
  row.names = FALSE
)

message("Generating matched single-cell IMC-Xenium score correlations")
imc_score_lookup <- tibble(
  KEY = colnames(spe7),
  imc_glycolysis = unname(glycolysis_score[colnames(spe7)]),
  imc_tca = unname(tca_score[colnames(spe7)])
)
cor_base <- energy %>%
  filter(imc_matched, !is.na(KEY), nzchar(KEY)) %>%
  inner_join(imc_score_lookup, by = "KEY")

fib_cor <- bind_rows(
  cor_base %>% transmute(
    cell,
    donor,
    metfiblabel,
    pathway = "Glycolysis",
    imc_score = imc_glycolysis,
    xenium_score = glycolysis_scmetabolism
  ),
  cor_base %>% transmute(
    cell,
    donor,
    metfiblabel,
    pathway = "TCA / OXPHOS",
    imc_score = imc_tca,
    xenium_score = tca_scmetabolism
  )
) %>%
  filter(metfiblabel %in% fib_met_levels) %>%
  mutate(metfiblabel = factor(metfiblabel, levels = fib_met_levels))
fib_cor_result <- make_correlation(
  fib_cor,
  "metfiblabel",
  fib_met_levels,
  fib_colors,
  "Fibroblast IMC protein vs Xenium scMetabolism scores"
)
save_plot_set(
  fib_cor_result$plot,
  "fibmet_imc_scmetabolism_correlation_seven_donors",
  7.2,
  4.8
)
write.csv(
  fib_cor,
  file.path(table_dir, "fibmet_imc_scmetabolism_correlation_cells_seven_donors.csv"),
  row.names = FALSE
)
write.csv(
  fib_cor_result$stats,
  file.path(table_dir, "fibmet_imc_scmetabolism_correlation_stats_seven_donors.csv"),
  row.names = FALSE
)

ec_cor <- bind_rows(
  cor_base %>% transmute(
    cell,
    donor,
    metEClabel,
    pathway = "Glycolysis",
    imc_score = imc_glycolysis,
    xenium_score = glycolysis_scmetabolism
  ),
  cor_base %>% transmute(
    cell,
    donor,
    metEClabel,
    pathway = "TCA / OXPHOS",
    imc_score = imc_tca,
    xenium_score = tca_scmetabolism
  )
) %>%
  filter(metEClabel %in% ec_met_levels) %>%
  mutate(metEClabel = factor(metEClabel, levels = ec_met_levels))
ec_cor_result <- make_correlation(
  ec_cor,
  "metEClabel",
  ec_met_levels,
  ec_met_colors,
  "Endothelial IMC protein vs Xenium scMetabolism scores"
)
save_plot_set(
  ec_cor_result$plot,
  "endmet_imc_scmetabolism_correlation_seven_donors",
  7.2,
  4.8
)
write.csv(
  ec_cor,
  file.path(table_dir, "endmet_imc_scmetabolism_correlation_cells_seven_donors.csv"),
  row.names = FALSE
)
write.csv(
  ec_cor_result$stats,
  file.path(table_dir, "endmet_imc_scmetabolism_correlation_stats_seven_donors.csv"),
  row.names = FALSE
)

panel_audit <- tibble(
  population = c(
    "Joined Xenium",
    "IMC FibMET",
    "Paired Xenium FibMET",
    "IMC EndMET",
    "Paired Xenium EndMET",
    "ExcludedValidation paired FibMET",
    "ExcludedValidation paired EndMET"
  ),
  n_cells = c(
    ncol(sfe7),
    nrow(fib_imc),
    length(fib_paired_cells),
    nrow(ec_imc),
    sum(as.character(sfe7$metEClabel) %in% ec_met_levels),
    sum(as.character(sfe7$SampleId) == "ExcludedValidation" & as.character(sfe7$metfiblabel) %in% fib_met_levels),
    sum(as.character(sfe7$SampleId) == "ExcludedValidation" & as.character(sfe7$metEClabel) %in% ec_met_levels)
  )
)
write.csv(
  panel_audit,
  file.path(panel_dir, "seven_donor_panel_cell_audit.csv"),
  row.names = FALSE
)

cat("\nSeven-donor panel cell audit\n")
print(panel_audit)
cat("\nSaved manuscript panels in:", normalizePath(panel_dir, winslash = "/"), "\n")
