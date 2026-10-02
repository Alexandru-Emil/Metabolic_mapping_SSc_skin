# Scientific workflow derived from Revision/Nan annotation/Data/RMD files/MintFlow/07_main_mintflow_manuscript_figure/redo_main_mintflow_cached_other_ec_labels_fast.R
# Run from the repository root with METABOLIC_INPUT_DIR set to the external project data.
source("data_preprocessing/common/config.R")
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

out_dir <- project_path("Revision/Nan annotation/Data/RMD files/MintFlow/07_main_mintflow_manuscript_figure")
full_script <- file.path(out_dir, "redo_main_mintflow_cached_other_ec_labels.R")
final_dir_candidates <- c(
  project_path("external/perturb_met_hi_fib_ec_with_frequency_dotplots"),
  project_path("Revision/Nan annotation/Data/RMD files/MintFlow/05_final_perturbation_report_and_paper_plots/perturb_met_hi_fib_ec_with_frequency_dotplots")
)
obs_candidates <- c(
  project_path("external/obs.csv"),
  project_path("Revision/Nan annotation/Data/RMD files/MintFlow/01_anndata_exports/anndata_bundle_combined_fib_ec/obs.csv")
)
cache_dir <- file.path(out_dir, ".plot_cache")
cache_version <- "other_ec_labels_d1d2_v1"
force_replot <- identical(Sys.getenv("MINTFLOW_FORCE_REPLOT"), "1")

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

final_dir <- first_existing_dir(final_dir_candidates, "ec_cell_level_perturbation_results.csv")
obs_path <- first_existing_file(obs_candidates)
cell_path <- file.path(final_dir, "ec_cell_level_perturbation_results.csv")
paired_path <- file.path(final_dir, "paper_replacement_ec_paired_probability_delta.csv")
plot_inputs <- c(cell_path, obs_path, paired_path[file.exists(paired_path)])

obs_df <- read_csv(obs_path, show_col_types = FALSE) %>%
  mutate(
    Patient = as.character(Patient),
    SampleId = as.character(SampleId),
    sample_id = as.character(sample_id),
    donor_id = case_when(
      !is.na(Patient) & Patient != "" & !is.na(SampleId) & SampleId != "" ~ paste0(Patient, " (", SampleId, ")"),
      !is.na(SampleId) & SampleId != "" ~ SampleId,
      !is.na(sample_id) & sample_id != "" ~ sample_id,
      TRUE ~ "unknown_sample"
    )
  )

donor_levels <- obs_df %>%
  distinct(donor_id) %>%
  arrange(donor_id) %>%
  pull(donor_id)

expected <- bind_rows(
  lapply(donor_levels, function(donor) {
    slug <- slugify(donor)
    tibble(
      stem = c(
        paste0("spatial_response_by_donor_", slug),
        paste0("two_layer_spatial_by_donor_", slug)
      ),
      extra = list(
        list(donor = donor, plot = "standalone_response"),
        list(donor = donor, plot = "two_layer")
      )
    )
  }),
  tibble(
    stem = c(
      "spatial_zoomed_met_hi_fib_niche_examples",
      "panel_A_perturbation_schematic",
      "panel_B_replacement_ec_state_frequency",
      "panel_C_per_cell_probability_shift_by_original_ec_state",
      "panel_D_baseline_niche_and_ec_response",
      "panel_D_spatial_delta_map_faceted_by_original_ec_state",
      "supplemental_two_layer_spatial_niche_and_response",
      "main_mintflow_replacement_4panel_figure"
    ),
    extra = list(
      list(plot = "zoom"),
      list(plot = "a"),
      list(plot = "b"),
      list(plot = "c"),
      list(plot = "d_main"),
      list(plot = "d_legacy_filename"),
      list(plot = "supp_two_layer"),
      list(plot = "combined")
    )
  )
)

all_current <- all(vapply(
  seq_len(nrow(expected)),
  function(i) cache_current(expected$stem[[i]], plot_inputs, expected$extra[[i]]),
  logical(1)
))

html_file <- file.path(out_dir, "main_mintflow_replacement_4panel_report.html")
if (all_current && file.exists(html_file)) {
  message("fast cache hit: all plot outputs are current")
  message(html_file)
  message(file.path(out_dir, "main_mintflow_replacement_4panel_figure.pdf"))
  message(file.path(out_dir, "panel_D_baseline_niche_and_ec_response.pdf"))
  quit(save = "no", status = 0)
}

message("cache miss: running full generator")
source(full_script, local = FALSE)
