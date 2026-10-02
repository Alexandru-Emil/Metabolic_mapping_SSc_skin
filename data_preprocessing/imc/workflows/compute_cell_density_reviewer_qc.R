# Scientific workflow derived from Revision/Nan annotation/Data/RMD files/reviewer_QC_analyses_20260729/02_IMC_cell_density_control/compute_cell_density_reviewer_qc.R
# Run from the repository root with METABOLIC_INPUT_DIR set to the external project data.
source("data_preprocessing/common/config.R")
options(stringsAsFactors = FALSE, width = 180)

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(dplyr)
  library(tidyr)
  library(FNN)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

set.seed(20260729)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x[1])) y else x
}

script_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
script_dir <- normalizePath(
  dirname(script_file %||% "imc_cell_density_reviewer_qc/compute_cell_density_reviewer_qc.R"),
  winslash = "/",
  mustWork = FALSE
)
output_dir <- file.path(script_dir, "outputs")
plot_dir <- file.path(output_dir, "plots")
table_dir <- file.path(output_dir, "tables")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

spe_path <- paste0(
  paste0(project_path(""), "/"),
  "IMC_metabolic_analysis_AEM/Datasets/spe20250202205602.rds"
)
stopifnot(file.exists(spe_path))

spe <- readRDS(spe_path)
required_fields <- c(
  "sample_id", "patient_id", "Donor", "Progression_skin", "celltype",
  "Glycolysis_z_c", "TCA_z_c", "width_px", "height_px"
)
missing_fields <- setdiff(required_fields, colnames(colData(spe)))
if (length(missing_fields)) {
  stop("Missing required colData fields: ", paste(missing_fields, collapse = ", "))
}
if (!"delaunay_80_new" %in% colPairNames(spe)) {
  stop("The finalized object does not contain the expected delaunay_80_new graph.")
}

cell_metadata <- as.data.frame(colData(spe))
coordinates <- spatialCoords(spe)
if (nrow(coordinates) != nrow(cell_metadata)) {
  stop("Spatial coordinate and cell metadata dimensions do not agree.")
}

graph <- colPair(spe, "delaunay_80_new")
edge_from <- S4Vectors::from(graph)
edge_to <- S4Vectors::to(graph)
if (any(cell_metadata$sample_id[edge_from] != cell_metadata$sample_id[edge_to])) {
  stop("The stored Delaunay graph contains cross-sample edges.")
}

edge_distance <- sqrt(rowSums(
  (coordinates[edge_from, , drop = FALSE] - coordinates[edge_to, , drop = FALSE])^2
))
graph_summary <- data.frame(
  graph = "delaunay_80_new",
  n_directed_edges = length(edge_from),
  n_cells = ncol(spe),
  mean_out_degree = length(edge_from) / ncol(spe),
  min_edge_distance_um = min(edge_distance),
  median_edge_distance_um = median(edge_distance),
  max_edge_distance_um = max(edge_distance),
  cross_sample_edges = 0L
)

delaunay_degree <- tabulate(edge_from, nbins = ncol(spe))

mean_distance_5nn <- rep(NA_real_, ncol(spe))
sample_indices <- split(seq_len(ncol(spe)), cell_metadata$sample_id)
for (indices in sample_indices) {
  if (length(indices) <= 5) {
    stop("A sample contains fewer than six cells; five-neighbor distance is undefined.")
  }
  current_knn <- FNN::get.knn(coordinates[indices, , drop = FALSE], k = 5)
  mean_distance_5nn[indices] <- rowMeans(current_knn$nn.dist)
}

cell_density_data <- bind_cols(
  cell_metadata,
  data.frame(
    delaunay_neighbors_80um = delaunay_degree,
    mean_distance_5nn_um = mean_distance_5nn
  )
)

met_hi_cells <- cell_density_data %>%
  filter(celltype == "stromal_2") %>%
  mutate(
    Progression_skin = factor(
      as.character(Progression_skin),
      levels = c("healthy", "progressive", "stable")
    )
  )

if (!nrow(met_hi_cells)) {
  stop("No stromal_2 / Met_hi_Fib cells were found.")
}

donor_density_scores <- met_hi_cells %>%
  group_by(Donor, patient_id, Progression_skin) %>%
  summarise(
    n_met_hi_fib = n(),
    mean_delaunay_neighbors_80um = mean(delaunay_neighbors_80um, na.rm = TRUE),
    mean_distance_5nn_um = mean(mean_distance_5nn_um, na.rm = TRUE),
    mean_glycolysis_score = mean(Glycolysis_z_c, na.rm = TRUE),
    mean_tca_oxphos_score = mean(TCA_z_c, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Progression_skin, Donor)

image_dimensions <- cell_density_data %>%
  distinct(Donor, patient_id, Progression_skin, sample_id, width_px, height_px)
if (any(table(image_dimensions$Donor) != 1L)) {
  stop("Each donor must have exactly one image-dimension record.")
}

donor_cell_counts <- cell_density_data %>%
  count(Donor, patient_id, Progression_skin, sample_id, name = "n_analyzed_cells") %>%
  left_join(
    image_dimensions,
    by = c("Donor", "patient_id", "Progression_skin", "sample_id")
  ) %>%
  mutate(
    image_area_mm2 = width_px * height_px / 1e6,
    analyzed_cells_per_mm2 = n_analyzed_cells / image_area_mm2,
    Progression_skin = factor(
      as.character(Progression_skin),
      levels = c("healthy", "progressive", "stable")
    )
  ) %>%
  arrange(Progression_skin, Donor)

bootstrap_spearman_ci <- function(x, y, iterations = 5000L, seed = 1L) {
  set.seed(seed)
  n <- length(x)
  boot <- replicate(iterations, {
    idx <- sample.int(n, n, replace = TRUE)
    if (
      length(unique(x[idx])) < 2 ||
      length(unique(y[idx])) < 2
    ) {
      return(NA_real_)
    }
    suppressWarnings(cor(x[idx], y[idx], method = "spearman"))
  })
  unname(quantile(boot, c(0.025, 0.975), na.rm = TRUE))
}

spearman_row <- function(data, density_col, score_col, min_met_hi_cells, scope) {
  current <- data %>%
    filter(n_met_hi_fib >= min_met_hi_cells)
  if (scope == "SSc only") {
    current <- current %>% filter(Progression_skin != "healthy")
  }
  test <- suppressWarnings(cor.test(
    current[[density_col]],
    current[[score_col]],
    method = "spearman",
    exact = FALSE
  ))
  ci <- bootstrap_spearman_ci(
    current[[density_col]],
    current[[score_col]],
    seed = 1000 + min_met_hi_cells +
      match(density_col, c(
        "mean_delaunay_neighbors_80um",
        "mean_distance_5nn_um"
      )) * 10 +
      match(score_col, c(
        "mean_glycolysis_score",
        "mean_tca_oxphos_score"
      )) * 100
  )
  data.frame(
    scope = scope,
    minimum_met_hi_fib_cells = min_met_hi_cells,
    density_metric = density_col,
    metabolic_score = score_col,
    n_donors = nrow(current),
    rho = unname(test$estimate),
    ci_low = ci[1],
    ci_high = ci[2],
    p_value = test$p.value
  )
}

density_columns <- c(
  "mean_delaunay_neighbors_80um",
  "mean_distance_5nn_um"
)
score_columns <- c(
  "mean_glycolysis_score",
  "mean_tca_oxphos_score"
)

primary_correlations <- bind_rows(lapply(density_columns, function(density_col) {
  bind_rows(lapply(score_columns, function(score_col) {
    spearman_row(
      donor_density_scores,
      density_col,
      score_col,
      min_met_hi_cells = 1L,
      scope = "All donors"
    )
  }))
})) %>%
  mutate(p_holm = p.adjust(p_value, method = "holm"))

threshold_correlations <- bind_rows(lapply(c(1L, 5L, 10L, 20L), function(minimum) {
  bind_rows(lapply(density_columns, function(density_col) {
    bind_rows(lapply(score_columns, function(score_col) {
      spearman_row(
        donor_density_scores,
        density_col,
        score_col,
        min_met_hi_cells = minimum,
        scope = "All donors"
      )
    }))
  }))
})) %>%
  group_by(minimum_met_hi_fib_cells) %>%
  mutate(p_holm = p.adjust(p_value, method = "holm")) %>%
  ungroup()

ssc_correlations <- bind_rows(lapply(density_columns, function(density_col) {
  bind_rows(lapply(score_columns, function(score_col) {
    spearman_row(
      donor_density_scores,
      density_col,
      score_col,
      min_met_hi_cells = 1L,
      scope = "SSc only"
    )
  }))
})) %>%
  mutate(p_holm = p.adjust(p_value, method = "holm"))

group_test <- function(data, value_col, analysis_name) {
  formula <- reformulate("Progression_skin", response = value_col)
  test <- kruskal.test(formula, data = data)
  data.frame(
    analysis = analysis_name,
    metric = value_col,
    n_healthy = sum(data$Progression_skin == "healthy"),
    n_progressive = sum(data$Progression_skin == "progressive"),
    n_stable = sum(data$Progression_skin == "stable"),
    statistic = unname(test$statistic),
    df = unname(test$parameter),
    p_value = test$p.value
  )
}

cell_count_group_tests <- bind_rows(
  group_test(donor_cell_counts, "n_analyzed_cells", "Analyzed cell count"),
  group_test(donor_cell_counts, "analyzed_cells_per_mm2", "Analyzed cells per mm2")
)

local_density_group_tests <- bind_rows(
  group_test(
    donor_density_scores,
    "mean_delaunay_neighbors_80um",
    "Mean Delaunay neighbor count around Met_hi_Fib"
  ),
  group_test(
    donor_density_scores,
    "mean_distance_5nn_um",
    "Mean distance to five nearest cells from Met_hi_Fib"
  )
)

group_medians <- bind_rows(
  donor_cell_counts %>%
    group_by(Progression_skin) %>%
    summarise(
      analysis = "Analyzed cell count",
      median = median(n_analyzed_cells),
      q25 = quantile(n_analyzed_cells, 0.25),
      q75 = quantile(n_analyzed_cells, 0.75),
      .groups = "drop"
    ),
  donor_cell_counts %>%
    group_by(Progression_skin) %>%
    summarise(
      analysis = "Analyzed cells per mm2",
      median = median(analyzed_cells_per_mm2),
      q25 = quantile(analyzed_cells_per_mm2, 0.25),
      q75 = quantile(analyzed_cells_per_mm2, 0.75),
      .groups = "drop"
    ),
  donor_density_scores %>%
    group_by(Progression_skin) %>%
    summarise(
      analysis = "Mean Delaunay neighbor count around Met_hi_Fib",
      median = median(mean_delaunay_neighbors_80um),
      q25 = quantile(mean_delaunay_neighbors_80um, 0.25),
      q75 = quantile(mean_delaunay_neighbors_80um, 0.75),
      .groups = "drop"
    ),
  donor_density_scores %>%
    group_by(Progression_skin) %>%
    summarise(
      analysis = "Mean distance to five nearest cells from Met_hi_Fib",
      median = median(mean_distance_5nn_um),
      q25 = quantile(mean_distance_5nn_um, 0.25),
      q75 = quantile(mean_distance_5nn_um, 0.75),
      .groups = "drop"
    )
)

write.csv(
  donor_density_scores,
  file.path(table_dir, "donor_met_hi_fib_density_and_scores.csv"),
  row.names = FALSE
)
write.csv(
  donor_cell_counts,
  file.path(table_dir, "donor_analyzed_cell_counts_and_density.csv"),
  row.names = FALSE
)
write.csv(
  primary_correlations,
  file.path(table_dir, "primary_donor_level_spearman_correlations.csv"),
  row.names = FALSE
)
write.csv(
  threshold_correlations,
  file.path(table_dir, "minimum_met_hi_cell_threshold_sensitivity.csv"),
  row.names = FALSE
)
write.csv(
  ssc_correlations,
  file.path(table_dir, "ssc_only_spearman_correlations.csv"),
  row.names = FALSE
)
write.csv(
  cell_count_group_tests,
  file.path(table_dir, "cell_count_group_kruskal_tests.csv"),
  row.names = FALSE
)
write.csv(
  local_density_group_tests,
  file.path(table_dir, "local_density_group_kruskal_tests.csv"),
  row.names = FALSE
)
write.csv(
  group_medians,
  file.path(table_dir, "group_medians_and_iqr.csv"),
  row.names = FALSE
)
write.csv(
  graph_summary,
  file.path(table_dir, "delaunay_graph_qc.csv"),
  row.names = FALSE
)

group_labels <- c(
  "healthy" = "Healthy",
  "progressive" = "Progressive SSc",
  "stable" = "Stable SSc"
)
group_colors <- c(
  "healthy" = "#6CC3F4",
  "progressive" = "#C03830",
  "stable" = "#F6BF93"
)
ink <- "#30363B"
grid_color <- "#E5E8EA"

theme_reviewer <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.18)),
      plot.subtitle = element_text(color = "#4B535A", size = rel(0.92)),
      plot.caption = element_text(color = "#596168", hjust = 0, size = rel(0.78)),
      strip.background = element_rect(fill = "#F3F4F5", color = "#AEB4B9"),
      strip.text = element_text(face = "bold", size = rel(0.93)),
      panel.grid.major = element_line(color = grid_color, linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = ink, linewidth = 0.4),
      axis.title = element_text(face = "bold"),
      legend.position = "top",
      legend.justification = "left",
      legend.box = "horizontal",
      plot.margin = margin(9, 15, 9, 9)
    )
}

theme_manuscript_density <- function(base_size = 8.4) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      strip.background = element_rect(
        fill = "#F3F4F5",
        color = "#AEB4B9",
        linewidth = 0.35
      ),
      strip.text = element_text(face = "bold", size = 8.2, lineheight = 0.95),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = grid_color, linewidth = 0.28),
      panel.border = element_rect(color = ink, linewidth = 0.45),
      axis.title = element_text(size = 8.6, face = "bold"),
      axis.text = element_text(size = 7.4, color = ink),
      axis.ticks = element_line(linewidth = 0.35),
      axis.ticks.length = grid::unit(1.3, "mm"),
      legend.position = "top",
      legend.justification = "center",
      legend.box = "horizontal",
      legend.box.just = "center",
      legend.title = element_text(size = 7.5),
      legend.text = element_text(size = 7.1),
      legend.key.height = grid::unit(3.4, "mm"),
      legend.key.width = grid::unit(4.0, "mm"),
      legend.spacing.x = grid::unit(1.3, "mm"),
      legend.margin = margin(0, 0, 2, 0),
      plot.margin = margin(4, 5, 4, 5, unit = "pt")
    )
}

save_plot <- function(plot, stem, width, height, dpi = 450) {
  ggsave(
    file.path(plot_dir, paste0(stem, ".pdf")),
    plot,
    width = width,
    height = height,
    device = cairo_pdf
  )
  ggsave(
    file.path(plot_dir, paste0(stem, ".svg")),
    plot,
    width = width,
    height = height,
    device = svglite::svglite
  )
  ggsave(
    file.path(plot_dir, paste0(stem, ".png")),
    plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
}

correlation_plot_data <- donor_density_scores %>%
  pivot_longer(
    cols = all_of(density_columns),
    names_to = "density_metric",
    values_to = "density_value"
  ) %>%
  pivot_longer(
    cols = all_of(score_columns),
    names_to = "metabolic_score",
    values_to = "score_value"
  ) %>%
  mutate(
    density_label = recode(
      density_metric,
      mean_delaunay_neighbors_80um = "Delaunay neighbors\n(maximum edge 80 µm)",
      mean_distance_5nn_um = "5-nearest-cell distance\n(µm)"
    ),
    score_label = recode(
      metabolic_score,
      mean_glycolysis_score = "Glycolysis",
      mean_tca_oxphos_score = "TCA/OXPHOS"
    ),
    panel_label = factor(
      paste(score_label, "vs", density_label),
      levels = c(
        "Glycolysis vs Delaunay neighbors\n(maximum edge 80 µm)",
        "Glycolysis vs 5-nearest-cell distance\n(µm)",
        "TCA/OXPHOS vs Delaunay neighbors\n(maximum edge 80 µm)",
        "TCA/OXPHOS vs 5-nearest-cell distance\n(µm)"
      )
    )
  )

correlation_annotations <- primary_correlations %>%
  mutate(
    density_label = recode(
      density_metric,
      mean_delaunay_neighbors_80um = "Delaunay neighbors\n(maximum edge 80 µm)",
      mean_distance_5nn_um = "5-nearest-cell distance\n(µm)"
    ),
    score_label = recode(
      metabolic_score,
      mean_glycolysis_score = "Glycolysis",
      mean_tca_oxphos_score = "TCA/OXPHOS"
    ),
    panel_label = factor(
      paste(score_label, "vs", density_label),
      levels = c(
        "Glycolysis vs Delaunay neighbors\n(maximum edge 80 µm)",
        "Glycolysis vs 5-nearest-cell distance\n(µm)",
        "TCA/OXPHOS vs Delaunay neighbors\n(maximum edge 80 µm)",
        "TCA/OXPHOS vs 5-nearest-cell distance\n(µm)"
      )
    ),
    annotation_x = if_else(
      density_metric == "mean_distance_5nn_um",
      Inf,
      -Inf
    ),
    annotation_y = case_when(
      density_metric == "mean_delaunay_neighbors_80um" ~ Inf,
      metabolic_score == "mean_glycolysis_score" ~ mean(range(
        donor_density_scores$mean_glycolysis_score,
        na.rm = TRUE
      )),
      TRUE ~ mean(range(
        donor_density_scores$mean_tca_oxphos_score,
        na.rm = TRUE
      ))
    ),
    annotation_hjust = if_else(
      density_metric == "mean_distance_5nn_um",
      1.05,
      -0.05
    ),
    annotation_vjust = if_else(
      density_metric == "mean_distance_5nn_um",
      0.5,
      1.08
    ),
    label = sprintf(
      "ρ = %.2f; p = %.3f\nn = %d donors",
      rho,
      p_value,
      n_donors
    )
  )

correlation_plot <- ggplot(
  correlation_plot_data,
  aes(x = density_value, y = score_value)
) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.48,
    linetype = 2,
    color = "#777F86"
  ) +
  geom_point(
    aes(fill = Progression_skin, size = n_met_hi_fib),
    shape = 21,
    color = ink,
    stroke = 0.38,
    alpha = 0.95
  ) +
  ggrepel::geom_text_repel(
    aes(label = Donor),
    size = 2.25,
    seed = 20260729,
    max.overlaps = Inf,
    min.segment.length = 0,
    box.padding = 0.16,
    point.padding = 0.11,
    segment.color = "#A3A9AE",
    segment.size = 0.20,
    force = 2.2,
    force_pull = 0.25,
    max.iter = 100000,
    max.time = 20
  ) +
  geom_label(
    data = correlation_annotations,
    aes(
      x = annotation_x,
      y = annotation_y,
      label = label,
      hjust = annotation_hjust,
      vjust = annotation_vjust
    ),
    inherit.aes = FALSE,
    size = 2.30,
    lineheight = 0.95,
    fill = "white",
    color = ink,
    linewidth = 0.18,
    label.padding = grid::unit(0.12, "lines"),
    label.r = grid::unit(0.08, "lines")
  ) +
  facet_wrap(vars(panel_label), ncol = 2, scales = "free") +
  scale_fill_manual(
    values = group_colors,
    breaks = names(group_labels),
    labels = unname(group_labels),
    name = NULL
  ) +
  scale_size_continuous(
    range = c(1.45, 3.30),
    trans = "sqrt",
    breaks = c(5, 20, 100, 500, 1000),
    name = "Met-hi Fib cells"
  ) +
  guides(
    size = guide_legend(
      order = 1,
      title.position = "left",
      title.theme = element_text(size = 7.5, vjust = 0.5)
    ),
    fill = guide_legend(
      order = 2,
      title.position = "left",
      override.aes = list(size = 2.2)
    )
  ) +
  labs(
    x = "Local-density metric",
    y = "Mean metabolic score"
  ) +
  theme_manuscript_density(8.4) +
  theme(
    strip.text = element_text(size = 8.0),
    panel.spacing = grid::unit(2.6, "mm"),
    legend.position = "top",
    legend.justification = "center",
    legend.box.just = "center",
    plot.margin = margin(3, 4, 3, 10, unit = "pt")
  )
save_plot(correlation_plot, "figure_1_density_metabolic_score_correlations", 7.1, 4.7)

count_plot_data <- donor_cell_counts %>%
  select(
    Donor,
    Progression_skin,
    n_analyzed_cells,
    analyzed_cells_per_mm2
  ) %>%
  pivot_longer(
    cols = c(n_analyzed_cells, analyzed_cells_per_mm2),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric_label = recode(
      metric,
      n_analyzed_cells = "Analyzed cells",
      analyzed_cells_per_mm2 = "Analyzed cells per mm²"
    )
  )

count_annotations <- cell_count_group_tests %>%
  mutate(
    metric_label = recode(
      metric,
      n_analyzed_cells = "Analyzed cells",
      analyzed_cells_per_mm2 = "Analyzed cells per mm²"
    ),
    label = sprintf("Kruskal-Wallis p = %.3f", p_value)
  )

cell_count_plot <- ggplot(
  count_plot_data,
  aes(x = Progression_skin, y = value, fill = Progression_skin)
) +
  geom_boxplot(
    width = 0.56,
    outlier.shape = NA,
    alpha = 0.42,
    color = ink,
    linewidth = 0.45
  ) +
  geom_point(
    position = position_jitter(width = 0.1, height = 0, seed = 20260729),
    shape = 21,
    size = 2.35,
    color = ink,
    stroke = 0.38
  ) +
  geom_label(
    data = count_annotations,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.04,
    vjust = 1.08,
    size = 2.45,
    fill = "white",
    color = ink,
    linewidth = 0.18,
    label.padding = grid::unit(0.12, "lines"),
    label.r = grid::unit(0.08, "lines")
  ) +
  facet_wrap(vars(metric_label), nrow = 1, scales = "free_y") +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_x_discrete(labels = c(
    healthy = "Healthy",
    progressive = "Progressive\nSSc",
    stable = "Stable\nSSc"
  )) +
  scale_y_continuous(
    labels = scales::label_comma(),
    expand = expansion(mult = c(0.05, 0.14))
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_manuscript_density(8.4) +
  theme(
    axis.text.x = element_text(
      size = 7.6,
      angle = 0,
      hjust = 0.5,
      vjust = 1,
      lineheight = 0.92
    ),
    strip.text = element_text(size = 8.3),
    panel.spacing = grid::unit(4.2, "mm"),
    plot.margin = margin(3, 4, 4, 5, unit = "pt")
  )
save_plot(cell_count_plot, "figure_2_analyzed_cell_count_and_global_density", 7.1, 2.8)

local_density_plot_data <- donor_density_scores %>%
  select(
    Donor,
    Progression_skin,
    n_met_hi_fib,
    mean_delaunay_neighbors_80um,
    mean_distance_5nn_um
  ) %>%
  pivot_longer(
    cols = c(mean_delaunay_neighbors_80um, mean_distance_5nn_um),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric_label = recode(
      metric,
      mean_delaunay_neighbors_80um = "Mean Delaunay neighbors\n(maximum edge 80 µm)",
      mean_distance_5nn_um = "Mean 5-nearest-cell distance\n(µm)"
    )
  )

local_density_annotations <- local_density_group_tests %>%
  mutate(
    metric_label = recode(
      metric,
      mean_delaunay_neighbors_80um = "Mean Delaunay neighbors\n(maximum edge 80 µm)",
      mean_distance_5nn_um = "Mean 5-nearest-cell distance\n(µm)"
    ),
    label = sprintf("Kruskal-Wallis p = %.3f", p_value)
  )

local_density_plot <- ggplot(
  local_density_plot_data,
  aes(x = Progression_skin, y = value, fill = Progression_skin)
) +
  geom_boxplot(
    width = 0.56,
    outlier.shape = NA,
    alpha = 0.35,
    color = ink,
    linewidth = 0.45
  ) +
  geom_point(
    aes(size = n_met_hi_fib),
    position = position_jitter(width = 0.1, height = 0, seed = 20260729),
    shape = 21,
    color = ink,
    stroke = 0.45
  ) +
  geom_label(
    data = local_density_annotations,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.1,
    size = 3.0,
    fill = "white",
    color = ink,
    linewidth = 0.2
  ) +
  facet_wrap(vars(metric_label), nrow = 1, scales = "free_y") +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_x_discrete(labels = group_labels) +
  scale_size_continuous(
    range = c(2.3, 4.8),
    trans = "sqrt",
    breaks = c(5, 20, 100, 500, 1000),
    name = "Met-hi Fib cells"
  ) +
  labs(
    title = "Local cellularity around Met-hi Fib by clinical group",
    subtitle = "Each point is one donor; point size indicates contributing Met-hi Fib cells",
    x = NULL,
    y = NULL
  ) +
  theme_reviewer(10.8) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    panel.spacing = unit(1.2, "lines")
  )
save_plot(local_density_plot, "figure_3_local_density_metrics_by_group", 10.2, 5.2)

threshold_plot_data <- threshold_correlations %>%
  mutate(
    minimum_label = factor(
      minimum_met_hi_fib_cells,
      levels = c(1, 5, 10, 20),
      labels = c("All donors", "≥5 cells", "≥10 cells", "≥20 cells")
    ),
    comparison = paste(
      recode(
        density_metric,
        mean_delaunay_neighbors_80um = "Delaunay neighbors",
        mean_distance_5nn_um = "5-NN distance"
      ),
      "vs",
      recode(
        metabolic_score,
        mean_glycolysis_score = "glycolysis",
        mean_tca_oxphos_score = "TCA/OXPHOS"
      )
    ),
    cell_label = sprintf("ρ %.2f\np %.2f", rho, p_value)
  )

threshold_plot <- ggplot(
  threshold_plot_data,
  aes(x = minimum_label, y = comparison, fill = rho)
) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = cell_label), size = 3.0, lineheight = 0.95, color = ink) +
  scale_fill_gradient2(
    low = "#3F78A8",
    mid = "white",
    high = "#D97706",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman ρ"
  ) +
  labs(
    title = "Sensitivity to minimum Met-hi Fib cell counts",
    subtitle = "Each column repeats the donor-level analysis after excluding sparsely represented donor means",
    x = "Minimum Met-hi Fib cells per donor",
    y = NULL,
    caption = "Raw two-sided Spearman p-values are displayed; Holm-adjusted values are available in the source table."
  ) +
  theme_reviewer(10.5) +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 25, hjust = 1)
  )
save_plot(threshold_plot, "figure_4_minimum_cell_threshold_sensitivity", 8.8, 5.7)

lower_plot <- (cell_count_plot | local_density_plot) +
  plot_layout(widths = c(0.95, 1.05))

combined_plot <- correlation_plot / lower_plot +
  plot_layout(heights = c(1.5, 1.0)) +
  plot_annotation(
    title = "Cell-density control analysis for IMC metabolic profiles",
    theme = theme(plot.title = element_text(face = "bold", size = 16))
  )
save_plot(combined_plot, "figure_5_combined_reviewer_density_control", 15.0, 13.0)

source_manifest <- data.frame(
  object_path = spe_path,
  object_class = class(spe)[1],
  n_cells = ncol(spe),
  n_donors = dplyr::n_distinct(cell_metadata$Donor),
  met_hi_fib_definition = "celltype == stromal_2",
  glycolysis_score = "Glycolysis_z_c",
  tca_oxphos_score = "TCA_z_c",
  neighbor_graph = "delaunay_80_new",
  nearest_neighbor_k = 5L
)
write.csv(source_manifest, file.path(table_dir, "source_manifest.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

primary_lines <- apply(primary_correlations, 1, function(row) {
  sprintf(
    "%s vs %s: rho = %.3f, 95%% bootstrap CI %.3f to %.3f, p = %.4f, Holm p = %.4f, n = %d donors",
    row[["density_metric"]],
    row[["metabolic_score"]],
    as.numeric(row[["rho"]]),
    as.numeric(row[["ci_low"]]),
    as.numeric(row[["ci_high"]]),
    as.numeric(row[["p_value"]]),
    as.numeric(row[["p_holm"]]),
    as.integer(row[["n_donors"]])
  )
})
summary_lines <- c(
  "IMC cell-density reviewer analysis",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  primary_lines,
  "",
  sprintf(
    "Analyzed cell count across groups: Kruskal-Wallis p = %.4f.",
    cell_count_group_tests$p_value[
      cell_count_group_tests$metric == "n_analyzed_cells"
    ]
  ),
  sprintf(
    "Analyzed cells per mm2 across groups: Kruskal-Wallis p = %.4f.",
    cell_count_group_tests$p_value[
      cell_count_group_tests$metric == "analyzed_cells_per_mm2"
    ]
  )
)
writeLines(summary_lines, file.path(output_dir, "analysis_summary.txt"))

message("Wrote reviewer-density outputs to: ", output_dir)
