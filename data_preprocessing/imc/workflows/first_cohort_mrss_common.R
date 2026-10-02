# Scientific workflow derived from Revision/first_cohort_mrss_common.R
# Run from the repository root with METABOLIC_INPUT_DIR set to the external project data.
source("data_preprocessing/common/config.R")
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

spe_path <- project_path("IMC_metabolic_analysis_AEM/Datasets/spe20250127174306.rds")
clinical_path <- project_path("Revision/Clinical data_metabolic_AEM_therapy_localmRss.xlsx")

# Unicode modifier letters render as true superscripts in HTML, SVG, PDF, and PNG.
superscript_state <- c(
  hi = "\u02b0\u2071",
  int = "\u2071\u207f\u1d57",
  low = "\u02e1\u1d52\u02b7"
)

hi <- unname(superscript_state["hi"])
int <- unname(superscript_state["int"])
low <- unname(superscript_state["low"])

cluster_labels <- c(
  stromal_1 = paste0("TCA/OXPHOS", hi, "\nPGC-1\u03b1", hi, "_Fib"),
  stromal_2 = paste0("Met", hi, "_Fib"),
  stromal_3 = paste0("TCA/OXPHOS", int, "\nPGC-1\u03b1", hi, "_Fib"),
  stromal_4 = paste0("pmTOR", hi, "CD98", hi, "_Fib"),
  stromal_5 = paste0("Met", low, "_Fib"),
  stromal_6 = paste0("Glycolysis", int, "\nCS", hi, "OGDH", low, "_Fib"),
  stromal_7 = paste0("Glut", hi, "_Fib"),
  stromal_8 = paste0("Glycolysis", int, "\nPGC-1\u03b1", low, "FAO", low, "_Fib"),
  CD31_1 = paste0("Met", low, "_EC"),
  CD31_2 = paste0("TCA/OXPHOS", hi, "\nPGC-1\u03b1", hi, "_EC"),
  CD31_3 = paste0("Met", hi, "_EC"),
  CD31_4 = paste0("pmTOR", hi, "CD98", hi, "_EC"),
  CD31_5 = paste0("Glycolysis", hi, "\nTCA/OXPHOS", int, "PPP", hi, "_EC"),
  CD45_CD68_1 = paste0("Met", hi, "_Mf"),
  CD45_CD68_2 = paste0("Met", low, "_Mf"),
  CD45_CD68_3 = paste0("Glycolysis", low, "\nTCA/OXPHOS", hi, "_Mf"),
  CD45_CD68_4 = paste0("Glycolysis", hi, "\nTCA/OXPHOS", hi, "_Mf")
)

cluster_labels_single_line <- gsub("\n", " ", cluster_labels, fixed = TRUE)

# Plot-only rich-text labels. ggtext renders these as true superscripts with a
# consistent font and baseline in every graphics device.
cluster_labels_html <- c(
  stromal_1 = "TCA/OXPHOS<sup>hi</sup><br>PGC-1\u03b1<sup>hi</sup>_Fib",
  stromal_2 = "Met<sup>hi</sup>_Fib",
  stromal_3 = "TCA/OXPHOS<sup>int</sup><br>PGC-1\u03b1<sup>hi</sup>_Fib",
  stromal_4 = "pmTOR<sup>hi</sup> CD98<sup>hi</sup>_Fib",
  stromal_5 = "Met<sup>low</sup>_Fib",
  stromal_6 = "Glycolysis<sup>int</sup><br>CS<sup>hi</sup> OGDH<sup>low</sup>_Fib",
  stromal_7 = "Glut<sup>hi</sup>_Fib",
  stromal_8 = "Glycolysis<sup>int</sup><br>PGC-1\u03b1<sup>low</sup> FAO<sup>low</sup>_Fib",
  CD31_1 = "Met<sup>low</sup>_EC",
  CD31_2 = "TCA/OXPHOS<sup>hi</sup><br>PGC-1\u03b1<sup>hi</sup>_EC",
  CD31_3 = "Met<sup>hi</sup>_EC",
  CD31_4 = "pmTOR<sup>hi</sup> CD98<sup>hi</sup>_EC",
  CD31_5 = "Glycolysis<sup>hi</sup><br>TCA/OXPHOS<sup>int</sup> PPP<sup>hi</sup>_EC",
  CD45_CD68_1 = "Met<sup>hi</sup>_Mf",
  CD45_CD68_2 = "Met<sup>low</sup>_Mf",
  CD45_CD68_3 = "Glycolysis<sup>low</sup><br>TCA/OXPHOS<sup>hi</sup>_Mf",
  CD45_CD68_4 = "Glycolysis<sup>hi</sup><br>TCA/OXPHOS<sup>hi</sup>_Mf"
)

cluster_labels_html_single_line <- gsub("<br>", " ", cluster_labels_html, fixed = TRUE)

label_audit <- tibble(
  cluster_id = names(cluster_labels),
  rendered_label = unname(cluster_labels_single_line),
  expected_states = c(
    "TCA/OXPHOS: hi; PGC-1\u03b1: hi",
    "Met: hi",
    "TCA/OXPHOS: int; PGC-1\u03b1: hi",
    "pmTOR: hi; CD98: hi",
    "Met: low",
    "Glycolysis: int; CS: hi; OGDH: low",
    "Glut: hi",
    "Glycolysis: int; PGC-1\u03b1: low; FAO: low",
    "Met: low",
    "TCA/OXPHOS: hi; PGC-1\u03b1: hi",
    "Met: hi",
    "pmTOR: hi; CD98: hi",
    "Glycolysis: hi; TCA/OXPHOS: int; PPP: hi",
    "Met: hi",
    "Met: low",
    "Glycolysis: low; TCA/OXPHOS: hi",
    "Glycolysis: hi; TCA/OXPHOS: hi"
  )
)

stopifnot(
  length(cluster_labels) == 17,
  length(cluster_labels_html) == 17,
  identical(names(cluster_labels), names(cluster_labels_html)),
  all(nzchar(cluster_labels)),
  all(grepl("<sup>(hi|int|low)</sup>", cluster_labels_html)),
  all(vapply(superscript_state, function(x) any(grepl(x, cluster_labels, fixed = TRUE)), logical(1)))
)

parse_local_score <- function(x) {
  x <- as.character(x)
  out <- suppressWarnings(as.numeric(x))
  needs_parse <- is.na(out) & grepl("[0-9]", x)
  out[needs_parse] <- suppressWarnings(as.numeric(
    sub(".*?([0-9]+(?:\\.[0-9]+)?).*", "\\1", x[needs_parse], perl = TRUE)
  ))
  out
}

exact_binary_permutation_p <- function(score, response, observed_rho) {
  n_high <- sum(score == max(score))
  allocations <- combn(seq_along(score), n_high)
  permuted <- apply(allocations, 2, function(high_indices) {
    permuted_score <- rep(min(score), length(score))
    permuted_score[high_indices] <- max(score)
    suppressWarnings(cor(permuted_score, response, method = "spearman"))
  })
  mean(
    abs(permuted) >= abs(observed_rho) - sqrt(.Machine$double.eps),
    na.rm = TRUE
  )
}

load_first_cohort_mrss_data <- function() {
  clinical_raw <- suppressWarnings(read_xlsx(clinical_path, sheet = "Summary"))
  clinical <- clinical_raw %>%
    transmute(
      patient_id = `Label in metabolic panel`,
      subtype = `Subtype (dc/lc)`,
      current_mRSS = suppressWarnings(as.numeric(`current mRSS`)),
      local_mRSS = parse_local_score(`Local skin score`)
    ) %>%
    filter(subtype != "healthy") %>%
    mutate(included_local_mRSS = !is.na(local_mRSS))

  spe <- readRDS(spe_path)
  cell_data <- as.data.frame(colData(spe)) %>%
    transmute(
      patient_id = patient_id,
      indication = indication,
      celltype = new_celllabels,
      metabolic_cluster = pheno_clusters_new2
    )

  list(
    clinical = clinical,
    cell_data = cell_data,
    neighbor_data = as.data.frame(spe$hasNeigh_df)
  )
}

summarise_correlations <- function(data, score_variable, group_variables, exact_binary = FALSE) {
  data %>%
    filter(!is.na(.data[[score_variable]])) %>%
    group_by(across(all_of(group_variables))) %>%
    group_modify(~{
      score <- .x[[score_variable]]
      response <- .x$frequency_percent
      test <- suppressWarnings(cor.test(score, response, method = "spearman", exact = FALSE))
      rho <- unname(test$estimate)
      tibble(
        n_donors = nrow(.x),
        rho = rho,
        p_spearman_raw = test$p.value,
        p_exact_permutation = if (exact_binary) {
          exact_binary_permutation_p(score, response, rho)
        } else {
          NA_real_
        },
        minimum_frequency_percent = min(response),
        maximum_frequency_percent = max(response)
      )
    }) %>%
    ungroup()
}

format_p <- function(p) {
  ifelse(p < 0.001, format(p, scientific = TRUE, digits = 2), sprintf("%.3f", p))
}

save_vector_and_png <- function(plot, stem, width, height, dpi = 400) {
  ggsave(paste0(stem, ".svg"), plot, width = width, height = height, units = "in", bg = "white")
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, units = "in", device = cairo_pdf, bg = "white")
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, units = "in", dpi = dpi, bg = "white")
}
