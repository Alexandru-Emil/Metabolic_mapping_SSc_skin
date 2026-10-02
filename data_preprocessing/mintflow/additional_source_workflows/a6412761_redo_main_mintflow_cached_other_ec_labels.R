source("data_preprocessing/common/config.R")
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(patchwork)
})

out_dir <- project_path("Revision/Nan annotation/Data/RMD files/MintFlow/07_main_mintflow_manuscript_figure")
final_dir_candidates <- c(
  project_path("external/perturb_met_hi_fib_ec_with_frequency_dotplots"),
  project_path("Revision/Nan annotation/Data/RMD files/MintFlow/05_final_perturbation_report_and_paper_plots/perturb_met_hi_fib_ec_with_frequency_dotplots")
)
obs_candidates <- c(
  project_path("external/obs.csv"),
  project_path("Revision/Nan annotation/Data/RMD files/MintFlow/01_anndata_exports/anndata_bundle_combined_fib_ec/obs.csv")
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cache_dir <- file.path(out_dir, ".plot_cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
force_replot <- identical(Sys.getenv("MINTFLOW_FORCE_REPLOT"), "1")
cache_version <- "other_ec_labels_d1d2_v1"

first_existing_dir <- function(paths, required_file) {
  hits <- paths[file.exists(file.path(paths, required_file))]
  if (!length(hits)) stop("Could not find required file: ", required_file)
  hits[[1]]
}

first_existing_file <- function(paths) {
  hits <- paths[file.exists(paths)]
  if (!length(hits)) stop("Could not find any of: ", paste(paths, collapse = "; "))
  hits[[1]]
}

slugify <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

input_signature <- function(paths, extra = list()) {
  paths <- normalizePath(paths[file.exists(paths)], mustWork = FALSE)
  data.frame(
    path = paths,
    mtime = as.numeric(file.info(paths)$mtime),
    size = file.info(paths)$size,
    stringsAsFactors = FALSE
  ) |>
    list(inputs = _, extra = extra, version = cache_version)
}

cache_current <- function(stem, paths, extra = list()) {
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  png_path <- file.path(out_dir, paste0(stem, ".png"))
  meta_path <- file.path(cache_dir, paste0(stem, ".rds"))
  if (force_replot || !file.exists(pdf_path) || !file.exists(png_path) || !file.exists(meta_path)) {
    return(FALSE)
  }
  identical(readRDS(meta_path), input_signature(paths, extra))
}

theme_paper <- function(base_size = 9) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0),
      legend.key = element_blank(),
      legend.background = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "grey70"),
      strip.text = element_text(face = "bold")
    )
}

save_plot_cached <- function(plot, stem, width, height, dpi = 350, inputs, extra = list()) {
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  png_path <- file.path(out_dir, paste0(stem, ".png"))
  meta_path <- file.path(cache_dir, paste0(stem, ".rds"))
  if (cache_current(stem, inputs, extra)) {
    message("cache hit: ", stem)
    return(invisible(c(pdf_path, png_path)))
  }
  ggsave(pdf_path, plot, width = width, height = height, units = "in", useDingbats = FALSE)
  ggsave(png_path, plot, width = width, height = height, units = "in", dpi = dpi)
  saveRDS(input_signature(inputs, extra), meta_path)
  message("wrote: ", stem)
  invisible(c(pdf_path, png_path))
}

final_dir <- first_existing_dir(final_dir_candidates, "ec_cell_level_perturbation_results.csv")
obs_path <- first_existing_file(obs_candidates)
cell_path <- file.path(final_dir, "ec_cell_level_perturbation_results.csv")
paired_path <- file.path(final_dir, "paper_replacement_ec_paired_probability_delta.csv")
plot_inputs <- c(cell_path, obs_path, paired_path[file.exists(paired_path)])

state_cols <- c("Met_hi_EC-like" = "#b95b14", "Other_EC-like" = "#2f8f47")
label_cols <- c("Met_hi_EC" = "#b95b14", "Other_EC" = "#2f8f47")
baseline_cols <- c(
  "Met_hi_Fib" = "#4a0b0b",
  "Other_Fib" = "grey76",
  "Met_hi_EC" = "#b95b14",
  "Other_EC" = "#2f8f47",
  "Other" = "grey88"
)

cell_df <- read_csv(cell_path, show_col_types = FALSE) %>%
  filter(scenario %in% c("baseline", "replace_Met_hi_Fib_with_Other_Fib")) %>%
  mutate(
    scenario_label = recode(
      scenario,
      baseline = "Baseline",
      replace_Met_hi_Fib_with_Other_Fib = "Replace Met_hi_Fib\nwith Other_Fib"
    ),
    scenario_label = factor(scenario_label, levels = c("Baseline", "Replace Met_hi_Fib\nwith Other_Fib")),
    predicted_ec_state_like = recode(
      predicted_ec_state,
      Met_hi_EC = "Met_hi_EC-like",
      Met_int_EC = "Other_EC-like"
    ),
    predicted_ec_state_like = factor(predicted_ec_state_like, levels = c("Met_hi_EC-like", "Other_EC-like"))
  )

if (file.exists(paired_path)) {
  paired_df <- read_csv(paired_path, show_col_types = FALSE)
} else {
  paired_df <- cell_df %>%
    select(cell_id, metEClabel, scenario, prob_Met_hi_EC) %>%
    pivot_wider(names_from = scenario, values_from = prob_Met_hi_EC) %>%
    rename(prob_baseline = baseline, prob_replace = replace_Met_hi_Fib_with_Other_Fib) %>%
    mutate(delta_prob = prob_replace - prob_baseline)
}

obs_df <- read_csv(obs_path, show_col_types = FALSE) %>%
  mutate(
    cell_id = as.character(cell_id),
    Patient = as.character(Patient),
    SampleId = as.character(SampleId),
    sample_id = as.character(sample_id),
    combined_met_label = as.character(combined_met_label),
    x_centroid = as.numeric(x_centroid),
    y_centroid = as.numeric(y_centroid),
    donor_id = case_when(
      !is.na(Patient) & Patient != "" & !is.na(SampleId) & SampleId != "" ~ paste0(Patient, " (", SampleId, ")"),
      !is.na(SampleId) & SampleId != "" ~ SampleId,
      !is.na(sample_id) & sample_id != "" ~ sample_id,
      TRUE ~ "unknown_sample"
    )
  )

if (!all(c("x_centroid", "y_centroid") %in% names(paired_df))) {
  paired_df <- paired_df %>%
    left_join(obs_df %>% select(cell_id, x_centroid, y_centroid), by = "cell_id")
} else if (any(is.na(paired_df$x_centroid)) || any(is.na(paired_df$y_centroid))) {
  paired_df <- paired_df %>%
    select(-any_of(c("x_centroid_obs", "y_centroid_obs"))) %>%
    left_join(obs_df %>% select(cell_id, x_centroid_obs = x_centroid, y_centroid_obs = y_centroid), by = "cell_id") %>%
    mutate(
      x_centroid = coalesce(x_centroid, x_centroid_obs),
      y_centroid = coalesce(y_centroid, y_centroid_obs)
    ) %>%
    select(-x_centroid_obs, -y_centroid_obs)
}

paired_df <- paired_df %>%
  select(-any_of(c("donor_id", "Patient", "SampleId", "sample_id", "combined_met_label"))) %>%
  left_join(
    obs_df %>% select(cell_id, donor_id, Patient, SampleId, sample_id, combined_met_label),
    by = "cell_id"
  ) %>%
  mutate(
    metEClabel_display = recode(as.character(metEClabel), Met_hi_EC = "Met_hi_EC", Met_int_EC = "Other_EC"),
    metEClabel_display = factor(metEClabel_display, levels = c("Met_hi_EC", "Other_EC"))
  )

write_csv(paired_df, file.path(out_dir, "main_figure_paired_ec_probability_delta.csv"))

panel_a <- ggplot() +
  annotate("rect", xmin = 0.05, xmax = 0.95, ymin = 0.76, ymax = 0.91, fill = "#f8ebe5", color = "#8a3d2a", linewidth = 0.45) +
  annotate("rect", xmin = 0.05, xmax = 0.95, ymin = 0.53, ymax = 0.68, fill = "#e8f1f8", color = "#315a7d", linewidth = 0.45) +
  annotate("rect", xmin = 0.05, xmax = 0.95, ymin = 0.30, ymax = 0.45, fill = "#eef4e8", color = "#4d6f39", linewidth = 0.45) +
  annotate("rect", xmin = 0.05, xmax = 0.95, ymin = 0.07, ymax = 0.22, fill = "#f5f0df", color = "#7f6332", linewidth = 0.45) +
  annotate("text", x = 0.5, y = 0.835, label = "IMC-defined Met_hi_Fib\nprotein-derived label", size = 3.1) +
  annotate("text", x = 0.5, y = 0.605, label = "Matched Xenium / MintFlow model\nRNA generative space", size = 3.1) +
  annotate("text", x = 0.5, y = 0.375, label = "Replace Met_hi_Fib\nwith Other_Fib", size = 3.1) +
  annotate("text", x = 0.5, y = 0.145, label = "Predict generated EC response\nMet_hi_EC-like probability", size = 3.1) +
  annotate("segment", x = 0.5, xend = 0.5, y = 0.76, yend = 0.69, arrow = arrow(length = grid::unit(0.12, "in")), linewidth = 0.5) +
  annotate("segment", x = 0.5, xend = 0.5, y = 0.53, yend = 0.46, arrow = arrow(length = grid::unit(0.12, "in")), linewidth = 0.5) +
  annotate("segment", x = 0.5, xend = 0.5, y = 0.30, yend = 0.23, arrow = arrow(length = grid::unit(0.12, "in")), linewidth = 0.5) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  labs(title = "A  Perturbation design") +
  theme_void(base_size = 9) +
  theme(plot.title = element_text(face = "bold", hjust = 0))

freq_df <- cell_df %>%
  count(scenario_label, predicted_ec_state_like, name = "n") %>%
  group_by(scenario_label) %>%
  mutate(total = sum(n), frequency = n / total) %>%
  ungroup()
write_csv(freq_df, file.path(out_dir, "main_figure_replacement_ec_state_frequency.csv"))

panel_b <- ggplot(freq_df, aes(x = scenario_label, y = frequency, fill = predicted_ec_state_like)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.68) +
  geom_text(
    aes(label = sprintf("%.1f%%", 100 * frequency)),
    position = position_dodge(width = 0.72),
    vjust = -0.35,
    size = 2.6
  ) +
  scale_fill_manual(values = state_cols, name = "Generated EC identity") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.74), expand = expansion(mult = c(0, 0.04))) +
  labs(
    title = "B  Replacement of Met_hi_Fib reduces Met_hi_EC-like generated ECs",
    x = NULL,
    y = "Predicted EC-state frequency"
  ) +
  theme_paper(9) +
  theme(legend.position = "right")

panel_c <- ggplot(paired_df, aes(x = metEClabel_display, y = delta_prob, fill = metEClabel_display)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(width = 0.16, size = 0.75, alpha = 0.42, color = "black") +
  scale_fill_manual(values = label_cols, guide = "none") +
  coord_cartesian(ylim = c(-1.05, 1.05)) +
  labs(
    title = "C  Per-cell EC probability shift",
    x = "Original IMC EC label",
    y = "Delta probability\nreplacement - baseline"
  ) +
  theme_paper(9)

all_cells_df <- obs_df %>%
  filter(!is.na(x_centroid), !is.na(y_centroid)) %>%
  mutate(
    baseline_label = case_when(
      combined_met_label == "Met_hi_Fib" ~ "Met_hi_Fib",
      combined_met_label == "Other_Fib" ~ "Other_Fib",
      combined_met_label == "Met_hi_EC" ~ "Met_hi_EC",
      combined_met_label == "Met_int_EC" ~ "Other_EC",
      TRUE ~ "Other"
    ),
    baseline_label = factor(
      baseline_label,
      levels = c("Met_hi_Fib", "Other_Fib", "Met_hi_EC", "Other_EC", "Other")
    )
  )

vmax <- max(0.1, as.numeric(quantile(abs(paired_df$delta_prob), probs = 0.99, na.rm = TRUE)))

make_two_layer_plot <- function(ec_df, tissue_df, title, point_scale = 1) {
  baseline_map <- ggplot(tissue_df, aes(x = x_centroid, y = y_centroid, color = baseline_label)) +
    geom_point(size = 0.72 * point_scale, alpha = 0.86) +
    scale_color_manual(values = baseline_cols, drop = FALSE, name = "Baseline label") +
    coord_equal() +
    scale_y_reverse() +
    labs(title = "D1  Baseline niche", x = "x", y = "y") +
    theme_paper(8) +
    theme(legend.position = "right")

  response_map <- ggplot() +
    geom_point(
      data = tissue_df %>% filter(combined_met_label != "Met_hi_Fib"),
      aes(x = x_centroid, y = y_centroid),
      color = "grey88",
      size = 0.32 * point_scale,
      alpha = 0.45
    ) +
    geom_point(
      data = tissue_df %>% filter(combined_met_label == "Met_hi_Fib"),
      aes(x = x_centroid, y = y_centroid),
      shape = 21,
      color = "grey10",
      fill = NA,
      size = 0.92 * point_scale,
      stroke = 0.34 * point_scale,
      alpha = 0.9
    ) +
    geom_point(
      data = ec_df,
      aes(x = x_centroid, y = y_centroid, color = delta_prob),
      size = 1.45 * point_scale,
      alpha = 0.98
    ) +
    scale_color_gradient2(
      low = "#2b59c3",
      mid = "grey94",
      high = "#c9332b",
      midpoint = 0,
      limits = c(-vmax, vmax),
      oob = scales::squish,
      name = "Delta Met_hi_EC-like\nprobability"
    ) +
    coord_equal() +
    scale_y_reverse() +
    labs(
      title = "D2  EC response after replacement",
      subtitle = "EC color: replacement - baseline. Black outlines: former Met_hi_Fib positions.",
      x = "x",
      y = "y"
    ) +
    theme_paper(8) +
    theme(
      legend.position = "right",
      plot.subtitle = element_text(size = 7.2, color = "grey25")
    )

  (baseline_map | response_map) +
    plot_annotation(
      title = title,
      theme = theme(plot.title = element_text(face = "bold", size = 10, hjust = 0))
    )
}

make_spatial_response_plot <- function(ec_df, tissue_df, title, point_scale = 1) {
  ggplot() +
    geom_point(
      data = tissue_df %>% filter(combined_met_label != "Met_hi_Fib"),
      aes(x = x_centroid, y = y_centroid),
      color = "grey88",
      size = 0.32 * point_scale,
      alpha = 0.55
    ) +
    geom_point(
      data = tissue_df %>% filter(combined_met_label == "Met_hi_Fib"),
      aes(x = x_centroid, y = y_centroid),
      shape = 21,
      size = 0.74 * point_scale,
      stroke = 0.22 * point_scale,
      color = "grey10",
      fill = NA,
      alpha = 0.86
    ) +
    geom_point(
      data = ec_df,
      aes(x = x_centroid, y = y_centroid, color = delta_prob),
      size = 1.45 * point_scale,
      alpha = 0.98
    ) +
    facet_wrap(~ metEClabel_display, nrow = 1, drop = FALSE, labeller = labeller(metEClabel_display = function(x) paste("Original", x))) +
    scale_color_gradient2(
      low = "#2b59c3",
      mid = "grey94",
      high = "#c9332b",
      midpoint = 0,
      limits = c(-vmax, vmax),
      oob = scales::squish,
      name = "EC color:\nDelta Met_hi_EC-like\nprobability"
    ) +
    coord_equal() +
    scale_y_reverse() +
    labs(
      title = title,
      subtitle = "Blue: lower after replacement; red: higher after replacement. Black outlines: former Met_hi_Fib. Light gray: other cells.",
      x = "x",
      y = "y"
    ) +
    theme_paper(8) +
    theme(
      legend.position = "right",
      plot.subtitle = element_text(size = 7.5, color = "grey25")
    )
}

panel_d <- make_two_layer_plot(
  ec_df = paired_df,
  tissue_df = all_cells_df,
  title = "D  Baseline niche and EC response to Met_hi_Fib replacement",
  point_scale = 1
)

donor_levels <- all_cells_df %>%
  distinct(donor_id) %>%
  arrange(donor_id) %>%
  pull(donor_id)

donor_outputs <- lapply(donor_levels, function(donor) {
  donor_slug <- slugify(donor)
  tissue_d <- all_cells_df %>% filter(donor_id == donor)
  ec_d <- paired_df %>% filter(donor_id == donor)

  response_stem <- paste0("spatial_response_by_donor_", donor_slug)
  response_plot <- make_spatial_response_plot(
    ec_df = ec_d,
    tissue_df = tissue_d,
    title = paste0("Spatial EC response to Met_hi_Fib replacement: ", donor),
    point_scale = 1.25
  )
  save_plot_cached(response_plot, response_stem, 11.4, 5.7, dpi = 500, inputs = plot_inputs, extra = list(donor = donor, plot = "standalone_response"))

  two_layer_stem <- paste0("two_layer_spatial_by_donor_", donor_slug)
  two_layer_plot <- make_two_layer_plot(
    ec_df = ec_d,
    tissue_df = tissue_d,
    title = paste0("Baseline niche and EC response: ", donor),
    point_scale = 1.25
  )
  save_plot_cached(two_layer_plot, two_layer_stem, 12.6, 5.6, dpi = 500, inputs = plot_inputs, extra = list(donor = donor, plot = "two_layer"))

  tibble(
    donor_id = donor,
    tab_id = paste0("tab_", donor_slug),
    response_png = paste0(response_stem, ".png"),
    response_pdf = paste0(response_stem, ".pdf"),
    two_layer_png = paste0(two_layer_stem, ".png"),
    two_layer_pdf = paste0(two_layer_stem, ".pdf"),
    n_cells = nrow(tissue_d),
    n_met_hi_fib = sum(tissue_d$combined_met_label == "Met_hi_Fib", na.rm = TRUE),
    n_ec = nrow(ec_d)
  )
}) %>%
  bind_rows()

write_csv(donor_outputs, file.path(out_dir, "spatial_map_outputs_by_donor.csv"))

zoom_html <- ""
zoom_centers <- paired_df %>%
  filter(!is.na(distance_to_nearest_baseline_Met_hi_Fib), !is.na(x_centroid), !is.na(y_centroid)) %>%
  arrange(distance_to_nearest_baseline_Met_hi_Fib, desc(abs(delta_prob))) %>%
  group_by(donor_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  arrange(distance_to_nearest_baseline_Met_hi_Fib) %>%
  slice_head(n = 3) %>%
  mutate(zoom_id = paste0("Zoom ", row_number(), ": ", donor_id))

if (nrow(zoom_centers) > 0) {
  zoom_half_width <- 120
  zoom_context <- lapply(seq_len(nrow(zoom_centers)), function(i) {
    center <- zoom_centers[i, ]
    all_cells_df %>%
      filter(
        donor_id == center$donor_id,
        x_centroid >= center$x_centroid - zoom_half_width,
        x_centroid <= center$x_centroid + zoom_half_width,
        y_centroid >= center$y_centroid - zoom_half_width,
        y_centroid <= center$y_centroid + zoom_half_width
      ) %>%
      mutate(zoom_id = center$zoom_id)
  }) %>% bind_rows()

  zoom_ec <- lapply(seq_len(nrow(zoom_centers)), function(i) {
    center <- zoom_centers[i, ]
    paired_df %>%
      filter(
        donor_id == center$donor_id,
        x_centroid >= center$x_centroid - zoom_half_width,
        x_centroid <= center$x_centroid + zoom_half_width,
        y_centroid >= center$y_centroid - zoom_half_width,
        y_centroid <= center$y_centroid + zoom_half_width
      ) %>%
      mutate(zoom_id = center$zoom_id)
  }) %>% bind_rows()

  zoom_plot <- ggplot() +
    geom_point(
      data = zoom_context %>% filter(combined_met_label != "Met_hi_Fib"),
      aes(x = x_centroid, y = y_centroid),
      color = "grey86",
      size = 0.8,
      alpha = 0.65
    ) +
    geom_point(
      data = zoom_context %>% filter(combined_met_label == "Met_hi_Fib"),
      aes(x = x_centroid, y = y_centroid),
      color = "grey8",
      fill = NA,
      shape = 21,
      size = 1.6,
      stroke = 0.25,
      alpha = 0.9
    ) +
    geom_point(
      data = zoom_ec,
      aes(x = x_centroid, y = y_centroid, color = delta_prob, shape = metEClabel_display),
      size = 2.4,
      alpha = 0.98
    ) +
    facet_wrap(~ zoom_id, scales = "free", nrow = 1) +
    scale_color_gradient2(
      low = "#2b59c3",
      mid = "grey94",
      high = "#c9332b",
      midpoint = 0,
      limits = c(-vmax, vmax),
      oob = scales::squish,
      name = "Delta Met_hi_EC-like\nprobability"
    ) +
    scale_shape_manual(values = c("Met_hi_EC" = 16, "Other_EC" = 17), name = "Original EC label") +
    scale_y_reverse() +
    labs(
      title = "Zoomed Met_hi_Fib niche examples",
      subtitle = "Windows are centered on ECs nearest to baseline Met_hi_Fib cells.",
      x = "x",
      y = "y"
    ) +
    theme_paper(8) +
    theme(legend.position = "right")
  save_plot_cached(zoom_plot, "spatial_zoomed_met_hi_fib_niche_examples", 12.4, 4.4, dpi = 500, inputs = plot_inputs, extra = list(plot = "zoom"))
  write_csv(zoom_centers, file.path(out_dir, "spatial_zoomed_met_hi_fib_niche_centers.csv"))
  zoom_html <- "<h2>Zoomed Niche Examples</h2><h3>spatial_zoomed_met_hi_fib_niche_examples.png</h3><img src='spatial_zoomed_met_hi_fib_niche_examples.png'>"
}

save_plot_cached(panel_a, "panel_A_perturbation_schematic", 4.6, 3.8, dpi = 350, inputs = plot_inputs, extra = list(plot = "a"))
save_plot_cached(panel_b, "panel_B_replacement_ec_state_frequency", 6.4, 4.2, dpi = 350, inputs = plot_inputs, extra = list(plot = "b"))
save_plot_cached(panel_c, "panel_C_per_cell_probability_shift_by_original_ec_state", 4.5, 4.1, dpi = 350, inputs = plot_inputs, extra = list(plot = "c"))
save_plot_cached(panel_d, "panel_D_baseline_niche_and_ec_response", 12.8, 5.7, dpi = 450, inputs = plot_inputs, extra = list(plot = "d_main"))
save_plot_cached(panel_d, "panel_D_spatial_delta_map_faceted_by_original_ec_state", 12.8, 5.7, dpi = 450, inputs = plot_inputs, extra = list(plot = "d_legacy_filename"))
save_plot_cached(panel_d, "supplemental_two_layer_spatial_niche_and_response", 12.8, 5.7, dpi = 450, inputs = plot_inputs, extra = list(plot = "supp_two_layer"))

combined <- (panel_a | panel_b) / (panel_c | panel_d) +
  plot_layout(widths = c(0.76, 1.52), heights = c(0.84, 1.22))
save_plot_cached(combined, "main_mintflow_replacement_4panel_figure", 15.2, 10.6, dpi = 450, inputs = plot_inputs, extra = list(plot = "combined"))

donor_tabs <- paste0(
  "<div class='tabs'>",
  paste0(
    "<button class='tablink' onclick=\"openTab(event, '", donor_outputs$tab_id, "')\">",
    html_escape(donor_outputs$donor_id),
    "</button>",
    collapse = ""
  ),
  "</div>",
  paste0(
    "<div id='", donor_outputs$tab_id, "' class='tabcontent'>",
    "<h3>", html_escape(donor_outputs$donor_id), "</h3>",
    "<p>Cells: ", donor_outputs$n_cells,
    "; baseline Met_hi_Fib: ", donor_outputs$n_met_hi_fib,
    "; ECs with paired MintFlow probabilities: ", donor_outputs$n_ec, ".</p>",
    "<p><a href='", donor_outputs$two_layer_pdf, "'>Baseline niche plus EC response PDF</a> | ",
    "<a href='", donor_outputs$response_pdf, "'>Standalone EC response PDF</a></p>",
    "<h4>Baseline niche and EC response side by side</h4>",
    "<img src='", donor_outputs$two_layer_png, "'>",
    "<h4>Supplemental standalone EC response with full tissue context</h4>",
    "<img src='", donor_outputs$response_png, "'>",
    "</div>",
    collapse = ""
  )
)

html_file <- file.path(out_dir, "main_mintflow_replacement_4panel_report.html")
html <- paste0(
  "<!doctype html><html><head><meta charset='utf-8'>",
  "<title>Main MintFlow Met_hi_Fib Replacement Figure</title>",
  "<style>body{font-family:Arial,sans-serif;margin:32px;color:#222;line-height:1.45} img{border:1px solid #ddd;margin-bottom:24px;max-width:100%} code{background:#f4f4f4;padding:1px 4px}.tabs{display:flex;flex-wrap:wrap;gap:6px;margin:16px 0}.tablink{border:1px solid #bbb;background:#f4f4f4;padding:8px 12px;cursor:pointer}.tablink.active{background:#222;color:#fff}.tabcontent{display:none;border-top:1px solid #ddd;padding-top:12px}</style>",
  "<script>function openTab(evt,id){var i,t,b;t=document.getElementsByClassName('tabcontent');for(i=0;i<t.length;i++){t[i].style.display='none';}b=document.getElementsByClassName('tablink');for(i=0;i<b.length;i++){b[i].className=b[i].className.replace(' active','');}document.getElementById(id).style.display='block';evt.currentTarget.className+=' active';}window.onload=function(){var b=document.getElementsByClassName('tablink');if(b.length){b[0].click();}}</script>",
  "</head><body>",
  "<h1>Main MintFlow Met_hi_Fib Replacement Figure</h1>",
  "<p>This report reorganizes the key MintFlow result from <code>mintflow_met_hi_fib_perturbation_ec_report_with_frequency_dotplots.html</code> into a manuscript-focused 4-panel layout. The former <code>Met_int_EC</code> label is displayed as <code>Other_EC</code>. Panel D is split into D1 baseline niche and D2 EC response after replacement: former Met_hi_Fib positions are overlaid as black outlines, all other cells are faint gray, and ECs are colored by replacement minus baseline Met_hi_EC-like probability.</p>",
  "<p><b>Source folder:</b> <code>", html_escape(final_dir), "</code></p>",
  "<p><b>Coordinate/label source:</b> <code>", html_escape(obs_path), "</code></p>",
  "<p><b>Cache:</b> plot cache is stored in <code>", html_escape(cache_dir), "</code>. Set <code>MINTFLOW_FORCE_REPLOT=1</code> before running the script to force all plots to be redrawn.</p>",
  "<h2>Combined Figure</h2>",
  "<h3>main_mintflow_replacement_4panel_figure.png</h3><img src='main_mintflow_replacement_4panel_figure.png'>",
  "<h2>Separate Panels</h2>",
  "<h3>panel_A_perturbation_schematic.png</h3><img src='panel_A_perturbation_schematic.png'>",
  "<h3>panel_B_replacement_ec_state_frequency.png</h3><img src='panel_B_replacement_ec_state_frequency.png'>",
  "<h3>panel_C_per_cell_probability_shift_by_original_ec_state.png</h3><img src='panel_C_per_cell_probability_shift_by_original_ec_state.png'>",
  "<h3>panel_D_baseline_niche_and_ec_response.png</h3><img src='panel_D_baseline_niche_and_ec_response.png'>",
  "<h2>Per-Donor Spatial Maps</h2>",
  "<p>Each tab shows one donor/sample separately, avoiding overlap from plotting all donors in one coordinate system.</p>",
  donor_tabs,
  "<h2>Supplemental Two-Layer Spatial Map</h2>",
  "<p>This plot first shows the baseline tissue niche and then the EC response after Met_hi_Fib replacement.</p>",
  "<h3>supplemental_two_layer_spatial_niche_and_response.png</h3><img src='supplemental_two_layer_spatial_niche_and_response.png'>",
  zoom_html,
  "</body></html>"
)
writeLines(html, html_file)

message(html_file)
message(file.path(out_dir, "main_mintflow_replacement_4panel_figure.pdf"))
message(file.path(out_dir, "panel_D_baseline_niche_and_ec_response.pdf"))
message(file.path(out_dir, "spatial_map_outputs_by_donor.csv"))
