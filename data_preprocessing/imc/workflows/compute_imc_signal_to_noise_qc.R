# Scientific workflow derived from Revision/Nan annotation/Data/RMD files/reviewer_QC_analyses_20260729/01_IMC_signal_to_noise_QC/compute_imc_signal_to_noise_qc.R
# Run from the repository root with METABOLIC_INPUT_DIR set to the external project data.
source("data_preprocessing/common/config.R")
options(stringsAsFactors = FALSE, width = 180)

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(mclust)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

set.seed(220224)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x[1])) y else x
}

script_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
script_dir <- normalizePath(
  dirname(script_file %||% "imc_signal_to_noise_qc/compute_imc_signal_to_noise_qc.R"),
  winslash = "/",
  mustWork = FALSE
)
output_dir <- file.path(script_dir, "outputs")
plot_dir <- file.path(output_dir, "plots")
table_dir <- file.path(output_dir, "tables")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

first_path <- paste0(
  paste0(project_path(""), "/"),
  "IMC_metabolic_analysis_AEM/Datasets/spe20250202205602.rds"
)
second_path <- paste0(
  paste0(project_path("Revision/"), "/"),
  "Nan annotation/Data/IMC spe objects/New IMC/spe2025-12-17 umap.RData"
)

stopifnot(file.exists(first_path), file.exists(second_path))

spe_first <- readRDS(first_path)
second_env <- new.env(parent = emptyenv())
load(second_path, envir = second_env)
stopifnot(exists("spe", envir = second_env, inherits = FALSE))
spe_second <- second_env$spe

technical_markers <- c(
  "80ArAr", "127I", "131Xe", "134Xe", "138Ba", "190BCKG",
  "191Ir", "193Ir", "DNA1", "DNA2", "ICSK1", "ICSK2",
  "208Pb", "209Bi"
)

canonical_marker <- function(x) {
  replacements <- c(
    "Glut1" = "GLUT1",
    "NRF2_p" = "pNRF2",
    "TFAM" = "mtTFA",
    "PDGFRA" = "PDGFRa",
    "S6_p" = "pS6",
    "pRibosom_S6" = "pS6",
    "PGC-1a" = "PGC1a",
    "E-cadherin" = "E-cadherin",
    "E_cadherin" = "E-cadherin",
    "Collagen" = "COL1A1",
    "PFKL1" = "PFKL/PFKM",
    "PFKL1_PFKM" = "PFKL/PFKM",
    "Ki-67" = "Ki67"
  )
  out <- x
  matched <- match(x, names(replacements))
  out[!is.na(matched)] <- unname(replacements[matched[!is.na(matched)]])
  out
}

get_sample_id <- function(spe) {
  cd <- as.data.frame(colData(spe))
  candidates <- c("sample_id", "patient_id", "roi")
  selected <- candidates[candidates %in% names(cd)][1]
  if (is.na(selected)) {
    stop("No sample, patient, or ROI identifier found in colData.")
  }
  as.character(cd[[selected]])
}

fit_marker_snr <- function(
    transformed,
    raw,
    sample_id,
    marker,
    cohort,
    min_cells_per_component = 20L) {
  keep <- is.finite(transformed) & is.finite(raw) & !is.na(sample_id)
  transformed <- transformed[keep]
  raw <- raw[keep]
  sample_id <- sample_id[keep]

  fit <- tryCatch(
    suppressWarnings(Mclust(transformed, G = 2, verbose = FALSE)),
    error = function(e) NULL
  )

  if (is.null(fit) || length(unique(fit$classification)) != 2) {
    global <- tibble(
      cohort = cohort,
      marker = marker,
      marker_canonical = canonical_marker(marker),
      n_cells = length(raw),
      model_name = NA_character_,
      signal_component = NA_integer_,
      background_component = NA_integer_,
      n_signal = NA_integer_,
      n_background = NA_integer_,
      fraction_signal = NA_real_,
      signal_mean = NA_real_,
      background_mean = NA_real_,
      snr = NA_real_,
      log2_signal = NA_real_,
      log2_snr = NA_real_,
      fit_status = "mixture_fit_failed"
    )
    return(list(global = global, sample = tibble()))
  }

  component <- as.integer(fit$classification)
  component_means <- tapply(raw, component, mean, na.rm = TRUE)
  signal_component <- as.integer(names(which.max(component_means)))
  background_component <- as.integer(names(which.min(component_means)))
  is_signal <- component == signal_component
  is_background <- component == background_component

  signal_mean <- mean(raw[is_signal], na.rm = TRUE)
  background_mean <- mean(raw[is_background], na.rm = TRUE)
  snr <- if (is.finite(background_mean) && background_mean > 0) {
    signal_mean / background_mean
  } else {
    NA_real_
  }

  global <- tibble(
    cohort = cohort,
    marker = marker,
    marker_canonical = canonical_marker(marker),
    n_cells = length(raw),
    model_name = fit$modelName %||% NA_character_,
    signal_component = signal_component,
    background_component = background_component,
    n_signal = sum(is_signal),
    n_background = sum(is_background),
    fraction_signal = mean(is_signal),
    signal_mean = signal_mean,
    background_mean = background_mean,
    snr = snr,
    log2_signal = ifelse(signal_mean > 0, log2(signal_mean), NA_real_),
    log2_snr = ifelse(snr > 0, log2(snr), NA_real_),
    fit_status = "ok"
  )

  sample_levels <- unique(sample_id)
  sample_stats <- lapply(sample_levels, function(current_sample) {
    idx <- sample_id == current_sample
    cur_signal <- raw[idx & is_signal]
    cur_background <- raw[idx & is_background]
    cur_signal_mean <- if (length(cur_signal)) mean(cur_signal, na.rm = TRUE) else NA_real_
    cur_background_mean <- if (length(cur_background)) mean(cur_background, na.rm = TRUE) else NA_real_
    cur_snr <- if (
      is.finite(cur_signal_mean) &&
      is.finite(cur_background_mean) &&
      cur_background_mean > 0
    ) {
      cur_signal_mean / cur_background_mean
    } else {
      NA_real_
    }

    tibble(
      cohort = cohort,
      sample_id = current_sample,
      marker = marker,
      marker_canonical = canonical_marker(marker),
      n_cells = sum(idx),
      n_signal = length(cur_signal),
      n_background = length(cur_background),
      fraction_signal = length(cur_signal) / sum(idx),
      signal_mean = cur_signal_mean,
      background_mean = cur_background_mean,
      snr = cur_snr,
      log2_signal = ifelse(cur_signal_mean > 0, log2(cur_signal_mean), NA_real_),
      log2_snr = ifelse(cur_snr > 0, log2(cur_snr), NA_real_),
      sufficient_components =
        length(cur_signal) >= min_cells_per_component &&
        length(cur_background) >= min_cells_per_component
    )
  })

  list(global = global, sample = bind_rows(sample_stats))
}

analyze_cohort <- function(spe, cohort) {
  stopifnot("counts" %in% assayNames(spe), "asinh" %in% assayNames(spe))
  markers <- setdiff(rownames(spe), technical_markers)
  sample_id <- get_sample_id(spe)
  counts_mat <- assay(spe, "counts")
  transformed_mat <- assay(spe, "asinh")

  results <- lapply(markers, function(marker) {
    fit_marker_snr(
      transformed = as.numeric(transformed_mat[marker, ]),
      raw = as.numeric(counts_mat[marker, ]),
      sample_id = sample_id,
      marker = marker,
      cohort = cohort
    )
  })

  list(
    global = bind_rows(lapply(results, `[[`, "global")),
    sample = bind_rows(lapply(results, `[[`, "sample"))
  )
}

message("Fitting first-cohort marker mixtures...")
first_results <- analyze_cohort(spe_first, "First cohort")
message("Fitting second-cohort marker mixtures...")
second_results <- analyze_cohort(spe_second, "Second cohort")

global_snr <- bind_rows(first_results$global, second_results$global)
sample_snr <- bind_rows(first_results$sample, second_results$sample)

sample_summary <- sample_snr %>%
  filter(sufficient_components, is.finite(log2_snr)) %>%
  group_by(cohort, marker, marker_canonical) %>%
  summarise(
    n_samples_total = n_distinct(sample_id),
    n_samples_valid = n(),
    median_log2_snr = median(log2_snr),
    q25_log2_snr = quantile(log2_snr, 0.25),
    q75_log2_snr = quantile(log2_snr, 0.75),
    min_log2_snr = min(log2_snr),
    max_log2_snr = max(log2_snr),
    median_signal = median(signal_mean),
    median_background = median(background_mean),
    .groups = "drop"
  )

is_displayed_marker <- function(cohort, marker_canonical) {
  !(cohort == "Second cohort" & marker_canonical == "TCF21")
}

# Retain all fitted channels in the QC tables, but omit the second-cohort-only
# TCF21 channel from presentation figures.
plot_global_snr <- global_snr %>%
  filter(is_displayed_marker(cohort, marker_canonical))

plot_sample_snr <- sample_snr %>%
  filter(is_displayed_marker(cohort, marker_canonical))

plot_sample_summary <- sample_summary %>%
  filter(is_displayed_marker(cohort, marker_canonical))

cohort_summary <- global_snr %>%
  group_by(cohort) %>%
  summarise(
    n_antibody_channels = n(),
    n_successful_mixture_fits = sum(fit_status == "ok"),
    median_snr = median(snr, na.rm = TRUE),
    q25_snr = quantile(snr, 0.25, na.rm = TRUE),
    q75_snr = quantile(snr, 0.75, na.rm = TRUE),
    min_snr = min(snr, na.rm = TRUE),
    max_snr = max(snr, na.rm = TRUE),
    median_signal = median(signal_mean, na.rm = TRUE),
    median_background = median(background_mean, na.rm = TRUE),
    .groups = "drop"
  )

shared_global <- global_snr %>%
  select(cohort, marker_canonical, marker, signal_mean, background_mean, snr, log2_snr) %>%
  pivot_wider(
    names_from = cohort,
    values_from = c(marker, signal_mean, background_mean, snr, log2_snr),
    names_sep = "__"
  ) %>%
  filter(
    !is.na(`log2_snr__First cohort`),
    !is.na(`log2_snr__Second cohort`)
  )

cross_cohort_spearman <- suppressWarnings(cor(
  shared_global$`log2_snr__First cohort`,
  shared_global$`log2_snr__Second cohort`,
  method = "spearman",
  use = "complete.obs"
))

write.csv(global_snr, file.path(table_dir, "global_marker_snr.csv"), row.names = FALSE)
write.csv(sample_snr, file.path(table_dir, "per_sample_marker_snr.csv"), row.names = FALSE)
write.csv(sample_summary, file.path(table_dir, "per_marker_sample_snr_summary.csv"), row.names = FALSE)
write.csv(cohort_summary, file.path(table_dir, "cohort_snr_summary.csv"), row.names = FALSE)
write.csv(shared_global, file.path(table_dir, "cross_cohort_shared_marker_snr.csv"), row.names = FALSE)

source_manifest <- tibble(
  cohort = c("First cohort", "Second cohort"),
  object_path = c(first_path, second_path),
  object_class = c(class(spe_first)[1], class(spe_second)[1]),
  n_cells = c(ncol(spe_first), ncol(spe_second)),
  n_samples = c(
    n_distinct(get_sample_id(spe_first)),
    n_distinct(get_sample_id(spe_second))
  ),
  n_antibody_channels = c(
    length(setdiff(rownames(spe_first), technical_markers)),
    length(setdiff(rownames(spe_second), technical_markers))
  ),
  transformed_assay = "asinh",
  untransformed_assay = "counts"
)
write.csv(source_manifest, file.path(table_dir, "source_manifest.csv"), row.names = FALSE)

palette <- c("First cohort" = "#246A9A", "Second cohort" = "#D97706")
neutral_dark <- "#33383D"
neutral_mid <- "#8C959D"
neutral_light <- "#E5E8EA"

theme_qc <- function(base_size = 11.5) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = base_size + 2,
        margin = margin(b = 4)
      ),
      plot.subtitle = element_text(
        color = "#4A5157",
        size = base_size - 0.2,
        margin = margin(b = 7)
      ),
      plot.caption = element_text(
        color = "#596168",
        hjust = 0,
        size = base_size - 1.2,
        margin = margin(t = 4)
      ),
      strip.background = element_rect(fill = "#F3F4F5", color = "#B8BEC3"),
      strip.text = element_text(face = "bold", size = base_size),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#E7E9EB", linewidth = 0.3),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.7),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 0.5),
      legend.position = "top",
      legend.justification = "left",
      legend.margin = margin(0, 0, 2, 0),
      plot.title.position = "plot",
      plot.margin = margin(6, 8, 6, 6)
    )
}

theme_manuscript_panel <- function(base_size = 8.2) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      strip.background = element_rect(
        fill = "#F3F4F5",
        color = "#AEB4B9",
        linewidth = 0.35
      ),
      strip.text = element_text(face = "bold", size = 8.2),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#E5E8EA", linewidth = 0.28),
      panel.border = element_rect(color = "#44484C", linewidth = 0.45),
      axis.title = element_text(size = 8.2),
      axis.text = element_text(size = 7.3),
      axis.ticks = element_line(linewidth = 0.35),
      axis.ticks.length = grid::unit(1.4, "mm"),
      legend.position = "none",
      plot.margin = margin(4, 4, 4, 5, unit = "pt")
    )
}

save_plot <- function(plot, stem, width, height, dpi = 400) {
  ggsave(file.path(plot_dir, paste0(stem, ".pdf")), plot, width = width, height = height)
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

reference_base <- ggplot(
  plot_global_snr %>%
    filter(fit_status == "ok", is.finite(log2_signal), is.finite(log2_snr)),
  aes(x = log2_signal, y = log2_snr)
) +
  geom_point(
    aes(fill = cohort),
    shape = 21,
    size = 3,
    stroke = 0.55,
    color = neutral_dark
  ) +
  scale_fill_manual(values = palette, guide = "none")

reference_plot <- reference_base +
  ggrepel::geom_text_repel(
    aes(label = marker_canonical),
    size = 3.1,
    seed = 220224,
    max.overlaps = Inf,
    min.segment.length = 0,
    box.padding = 0.28,
    point.padding = 0.18,
    segment.color = "#A3A9AE",
    segment.size = 0.28,
    force = 1.5,
    force_pull = 0.5,
    max.time = 3
  ) +
  facet_wrap(vars(cohort), ncol = 1) +
  labs(
    title = "Signal intensity and SNR",
    subtitle = "Mixtures fitted on asinh intensities; component means and SNR calculated on counts",
    x = expression(log[2] * " mean signal intensity"),
    y = expression(log[2] * " signal-to-noise ratio")
  ) +
  theme_qc(11.5)
save_plot(reference_plot, "figure_1_reference_style_signal_vs_snr", 7.1, 10.3)

marker_order <- plot_global_snr %>%
  group_by(marker_canonical) %>%
  summarise(order_value = median(log2_snr, na.rm = TRUE), .groups = "drop") %>%
  arrange(order_value) %>%
  pull(marker_canonical)

interval_data <- plot_sample_summary %>%
  mutate(marker_canonical = factor(marker_canonical, levels = marker_order))

interval_plot <- ggplot(
  interval_data,
  aes(
    x = median_log2_snr,
    y = marker_canonical,
    xmin = q25_log2_snr,
    xmax = q75_log2_snr,
    color = cohort
  )
) +
  geom_vline(xintercept = 0, color = neutral_mid, linewidth = 0.45) +
  geom_errorbarh(
    position = position_dodge(width = 0.55),
    height = 0,
    linewidth = 0.95
  ) +
  geom_point(
    position = position_dodge(width = 0.55),
    size = 2.8
  ) +
  scale_color_manual(values = palette) +
  labs(
    title = "Per-sample SNR consistency",
    subtitle = "Point: median; line: interquartile range across samples",
    x = expression(log[2] * " signal-to-noise ratio"),
    y = NULL,
    color = NULL
  ) +
  theme_qc(11.5) +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = "top"
  )
save_plot(interval_plot, "figure_2_per_sample_snr_intervals", 7.1, 9.4)

dumbbell_data <- plot_global_snr %>%
  filter(fit_status == "ok") %>%
  mutate(
    marker_canonical = factor(marker_canonical, levels = marker_order),
    log2_background = ifelse(background_mean > 0, log2(background_mean), NA_real_)
  )

dumbbell_plot <- ggplot(dumbbell_data, aes(y = marker_canonical)) +
  geom_segment(
    aes(
      x = log2_background,
      xend = log2_signal,
      yend = marker_canonical
    ),
    color = neutral_light,
    linewidth = 1.4,
    lineend = "round"
  ) +
  geom_point(aes(x = log2_background), shape = 21, fill = "white", color = neutral_mid, size = 2.8) +
  geom_point(aes(x = log2_signal, fill = cohort), shape = 21, color = neutral_dark, size = 3.2, stroke = 0.55) +
  facet_wrap(vars(cohort), nrow = 1) +
  scale_fill_manual(values = palette, guide = "none") +
  labs(
    title = "Empirical signal and background",
    subtitle = "Open: lower-intensity component; filled: higher-intensity component",
    x = expression(log[2] * " mean untransformed single-cell intensity"),
    y = NULL
  ) +
  theme_qc(11.5) +
  theme(axis.text.y = element_text(size = 10))
save_plot(dumbbell_plot, "figure_3_signal_background_dumbbell", 7.1, 9.4)

# Manuscript-ready Figure S2 panels are built at their final physical size.
# This avoids shrinking full-width report plots after assembly.
manuscript_signal_data <- plot_global_snr %>%
  filter(fit_status == "ok", is.finite(log2_signal), is.finite(log2_snr))

manuscript_signal_labels <- manuscript_signal_data

manuscript_signal_plot <- ggplot(
  manuscript_signal_data,
  aes(x = log2_signal, y = log2_snr)
) +
  geom_point(
    aes(fill = cohort),
    shape = 21,
    size = 1.55,
    stroke = 0.38,
    color = neutral_dark
  ) +
  scale_fill_manual(values = palette, guide = "none") +
  ggrepel::geom_text_repel(
    data = manuscript_signal_labels,
    aes(label = marker_canonical),
    size = 1.95,
    seed = 220224,
    max.overlaps = Inf,
    min.segment.length = 0,
    box.padding = 0.10,
    point.padding = 0.07,
    segment.color = "#9EA5AA",
    segment.size = 0.20,
    force = 2.8,
    force_pull = 0.18,
    max.iter = 100000,
    max.time = 20
  ) +
  facet_wrap(vars(cohort), ncol = 1) +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.08))) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.08))) +
  labs(
    x = expression(log[2] * " mean signal intensity"),
    y = expression(log[2] * " signal-to-noise ratio")
  ) +
  theme_manuscript_panel(8.2) +
  theme(
    axis.text = element_text(size = 7.4),
    strip.text = element_text(size = 8.4),
    panel.spacing.y = grid::unit(2.2, "mm"),
    plot.margin = margin(6, 3, 4, 7, unit = "pt")
  )

manuscript_background_plot <- ggplot(
  dumbbell_data,
  aes(y = marker_canonical)
) +
  geom_segment(
    aes(
      x = log2_background,
      xend = log2_signal,
      yend = marker_canonical
    ),
    color = "#B9BEC2",
    linewidth = 0.62,
    lineend = "round"
  ) +
  geom_point(
    aes(x = log2_background),
    shape = 21,
    fill = "white",
    color = "#70777D",
    size = 1.70,
    stroke = 0.38
  ) +
  geom_point(
    aes(x = log2_signal, fill = cohort),
    shape = 21,
    color = neutral_dark,
    size = 1.88,
    stroke = 0.38
  ) +
  facet_grid(. ~ cohort) +
  scale_fill_manual(values = palette, guide = "none") +
  scale_x_continuous(
    breaks = c(-6, 0, 6),
    expand = expansion(mult = c(0.04, 0.05))
  ) +
  labs(
    x = expression(log[2] * " mean component intensity"),
    y = NULL
  ) +
  theme_manuscript_panel(8.2) +
  theme(
    axis.text.x = element_text(size = 6.8),
    axis.text.y = element_text(size = 6.4),
    strip.text = element_text(size = 7.8),
    panel.spacing.x = grid::unit(2.2, "mm"),
    plot.margin = margin(6, 4, 4, 3, unit = "pt")
  )

figure_s2_combined <- wrap_plots(
  manuscript_signal_plot,
  manuscript_background_plot,
  nrow = 1,
  widths = c(1.30, 1.00)
) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 11),
      plot.tag.position = c(0, 1)
    )
  )

save_plot(
  manuscript_signal_plot,
  "figure_S2A_manuscript_signal_vs_snr",
  4.00,
  6.10
)
save_plot(
  manuscript_background_plot,
  "figure_S2B_manuscript_signal_background",
  3.10,
  6.10
)
save_plot(
  figure_s2_combined,
  "figure_S2_manuscript_combined",
  7.1,
  6.20
)

heatmap_data <- plot_sample_snr %>%
  filter(sufficient_components, is.finite(log2_snr)) %>%
  mutate(
    marker_canonical = factor(marker_canonical, levels = rev(marker_order)),
    sample_short = ifelse(
      grepl("^Slide_Xenium_", sample_id),
      sub(
        "^Slide_Xenium_([0-9]+)_[0-9.]+_([0-9]+)$",
        "\\1-\\2",
        sample_id
      ),
      sample_id
    )
  )

heat_limits <- range(heatmap_data$log2_snr, na.rm = TRUE)

make_heatmap <- function(data, cohort_name) {
  ggplot(
    data %>% filter(cohort == cohort_name),
    aes(x = sample_short, y = marker_canonical, fill = log2_snr)
  ) +
    geom_tile(color = "white", linewidth = 0.15) +
    scale_fill_gradient(
      low = "white",
      high = palette[[cohort_name]],
      limits = heat_limits,
      na.value = "#E8EAEC"
    ) +
    labs(
      title = paste(cohort_name, "sample-level SNR"),
      subtitle = "Common, unclipped color scale across cohorts",
      x = NULL,
      y = NULL,
      fill = expression(log[2] * " SNR")
    ) +
    theme_qc(11.5) +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8.5),
      axis.text.y = element_text(size = 8.9),
      panel.grid = element_blank()
    ) +
    guides(
      fill = guide_colorbar(
        barheight = grid::unit(35, "mm"),
        barwidth = grid::unit(5, "mm")
      )
    )
}

first_heatmap <- make_heatmap(heatmap_data, "First cohort")
second_heatmap <- make_heatmap(heatmap_data, "Second cohort")
save_plot(first_heatmap, "figure_4a_first_cohort_sample_snr_heatmap", 7.1, 8.0)
save_plot(second_heatmap, "figure_4b_second_cohort_sample_snr_heatmap", 7.1, 8.0)

concordance_plot <- ggplot(
  shared_global,
  aes(x = `log2_snr__First cohort`, y = `log2_snr__Second cohort`)
) +
  geom_abline(intercept = 0, slope = 1, color = neutral_mid, linetype = 2, linewidth = 0.5) +
  geom_point(shape = 21, fill = "#D7A12E", color = neutral_dark, size = 3, stroke = 0.55) +
  ggrepel::geom_text_repel(
    aes(label = marker_canonical),
    size = 2.9,
    seed = 220224,
    max.overlaps = Inf,
    min.segment.length = 0,
    box.padding = 0.25,
    point.padding = 0.16,
    segment.color = "#A3A9AE",
    segment.size = 0.28
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    hjust = -0.08,
    vjust = 1.1,
    label = sprintf("Spearman rho = %.2f; n = %d", cross_cohort_spearman, nrow(shared_global)),
    size = 3.2
  ) +
  coord_equal() +
  labs(
    title = "Cross-cohort SNR concordance",
    x = expression("First cohort " * log[2] * " SNR"),
    y = expression("Second cohort " * log[2] * " SNR")
  ) +
  theme_qc(11.5) +
  theme(legend.position = "none")
save_plot(concordance_plot, "figure_5_cross_cohort_snr_concordance", 7.1, 6.6)

reference_panel <- reference_base +
  ggrepel::geom_text_repel(
    aes(label = marker_canonical),
    size = 2.65,
    seed = 220224,
    max.overlaps = Inf,
    min.segment.length = 0,
    box.padding = 0.2,
    point.padding = 0.12,
    segment.color = "#A3A9AE",
    segment.size = 0.24
  ) +
  facet_wrap(vars(cohort), nrow = 1) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = expression(log[2] * " mean signal intensity"),
    y = expression(log[2] * " signal-to-noise ratio")
  ) +
  theme_qc(9.5) +
  theme(plot.margin = margin(3, 4, 3, 3))
interval_panel <- interval_plot +
  labs(title = NULL, subtitle = NULL) +
  theme(plot.margin = margin(3, 4, 3, 3))

combined_overview <- reference_panel / interval_panel +
  plot_layout(heights = c(0.85, 1.25)) +
  plot_annotation(
    title = "IMC final-panel signal-to-noise quality control",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.tag = element_text(face = "bold", size = 12)
    )
  )
save_plot(combined_overview, "figure_6_combined_qc_overview", 8.0, 11.2)

fit_diagnostics <- global_snr %>%
  transmute(
    cohort,
    marker,
    marker_canonical,
    fit_status,
    n_cells,
    n_signal,
    n_background,
    fraction_signal,
    signal_mean,
    background_mean,
    snr,
    log2_snr,
    component_size_flag =
      !is.na(n_signal) &
      !is.na(n_background) &
      pmin(n_signal, n_background) < 100
  )
write.csv(fit_diagnostics, file.path(table_dir, "mixture_fit_diagnostics.csv"), row.names = FALSE)

session_lines <- capture.output(sessionInfo())
writeLines(session_lines, file.path(output_dir, "sessionInfo.txt"))

summary_lines <- c(
  "IMC signal-to-noise QC summary",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  paste0(
    cohort_summary$cohort, ": ",
    cohort_summary$n_antibody_channels, " antibody channels; ",
    cohort_summary$n_successful_mixture_fits, " successful mixture fits; median SNR ",
    formatC(cohort_summary$median_snr, digits = 2, format = "f"),
    " (IQR ",
    formatC(cohort_summary$q25_snr, digits = 2, format = "f"), "-",
    formatC(cohort_summary$q75_snr, digits = 2, format = "f"), ")."
  ),
  "",
  sprintf(
    "Shared-marker cross-cohort Spearman rho: %.3f (n = %d).",
    cross_cohort_spearman,
    nrow(shared_global)
  ),
  "",
  paste(
    "Metric: a two-component Gaussian mixture was fitted to asinh-transformed",
    "single-cell mean intensities for each antibody and cohort. The component",
    "with the larger untransformed mean was designated signal and the lower",
    "component empirical background/noise. SNR = mean(signal counts) /",
    "mean(background counts)."
  )
)
writeLines(summary_lines, file.path(output_dir, "analysis_summary.txt"))

message("Wrote outputs to: ", output_dir)
