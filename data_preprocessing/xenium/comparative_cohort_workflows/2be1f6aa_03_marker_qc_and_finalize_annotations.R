suppressPackageStartupMessages({
  library(SpatialFeatureExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(svglite)
})

options(stringsAsFactors = FALSE)

project_root <- normalizePath(
  file.path(getwd(), "seven_donor_reanalysis"),
  winslash = "/",
  mustWork = TRUE
)
output_dir <- file.path(project_root, "outputs")
plot_dir <- file.path(output_dir, "annotation_qc_plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

annotated_path <- file.path(output_dir, "sfe_seven_donors_annotated_final.rds")
joint_candidate_path <- file.path(
  output_dir,
  "sfe_seven_donors_annotated_joint_candidate.rds"
)

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

major_markers <- list(
  epithelial = c("KRT5", "CDH1"),
  fibroblast = c("COL1A1", "PDGFRA"),
  pericyte = c("RGS5", "PDGFRB"),
  endothelial = c("VWF", "ENG"),
  muscle = c("MYL9", "DES"),
  telocyte = c("FOXL1", "PLN"),
  `Schwann cell` = c("MPZ", "SCN7A"),
  `myeloid cell` = c("MRC1", "CD1C"),
  lymphocyte = c("CD2", "MS4A1")
)
fib_markers <- list(
  papillary_Fib = c("COL18A1", "APCDD1"),
  COMP_Fib = c("COMP", "DPP4"),
  COL8A1_Fib = c("COL8A1", "SFRP4"),
  PI16_Fib = c("PI16", "MFAP5"),
  CXCL12_Fib = c("CXCL12", "APOE"),
  CCL19_Fib = c("CCL19", "CH25H"),
  NGFR_Fib = c("NGFR", "SCN7A"),
  COCH_Fib = c("COCH", "TNN")
)
ec_markers <- list(
  cycling_EC = c("MKI67", "TOP2A"),
  ACKR1_EC = c("ACKR1", "SELP"),
  EPC = c("KDR", "APLNR", "CD34"),
  ACTA2_EC = c("ACTA2"),
  HEY1_EC = c("HEY1", "GJA5"),
  lymphatic_EC = c("LYVE1", "PDPN")
)

strategy_labels <- c(
  joint_harmony = "Joint Harmony",
  reference_mapping = "Reference mapping"
)

save_plot_set <- function(plot, stem, width, height) {
  ggsave(
    file.path(plot_dir, paste0(stem, ".svg")),
    plot,
    width = width,
    height = height,
    units = "in",
    bg = "white",
    device = svglite::svglite
  )
  ggsave(
    file.path(plot_dir, paste0(stem, ".pdf")),
    plot,
    width = width,
    height = height,
    units = "in",
    bg = "white",
    device = cairo_pdf
  )
  ggsave(
    file.path(plot_dir, paste0(stem, ".png")),
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}

marker_availability <- function(marker_sets, genes, compartment) {
  tibble(
    compartment = compartment,
    annotation = rep(names(marker_sets), lengths(marker_sets)),
    gene = unlist(marker_sets, use.names = FALSE),
    available = unlist(marker_sets, use.names = FALSE) %in% genes
  )
}

make_dot_summary <- function(
    expression,
    counts,
    labels,
    marker_sets,
    strategy,
    level_order) {
  marker_order <- unique(unlist(marker_sets, use.names = FALSE))
  marker_order <- marker_order[marker_order %in% rownames(expression)]
  keep <- !is.na(labels) & labels %in% level_order
  labels <- factor(labels[keep], levels = level_order)
  expression <- expression[marker_order, keep, drop = FALSE]
  counts <- counts[marker_order, keep, drop = FALSE]

  result <- bind_rows(lapply(levels(droplevels(labels)), function(label) {
    idx <- labels == label
    tibble(
      annotation = label,
      gene = marker_order,
      mean_expression = Matrix::rowMeans(expression[, idx, drop = FALSE]),
      percent_expressed = 100 * Matrix::rowMeans(counts[, idx, drop = FALSE] > 0),
      n_cells = sum(idx)
    )
  }))

  result %>%
    group_by(gene) %>%
    mutate(
      scaled_mean_expression = if (n() > 1L && sd(mean_expression) > 0) {
        as.numeric(base::scale(mean_expression))
      } else {
        0
      }
    ) %>%
    ungroup() %>%
    mutate(
      strategy = strategy,
      strategy_label = unname(strategy_labels[strategy]),
      annotation = factor(annotation, levels = rev(level_order)),
      gene = factor(gene, levels = marker_order)
    )
}

signature_specificity <- function(
    expression,
    labels,
    marker_sets,
    strategy,
    compartment) {
  labels <- as.character(labels)
  marker_sets <- lapply(marker_sets, function(x) {
    x[x %in% rownames(expression)]
  })
  marker_sets <- marker_sets[lengths(marker_sets) > 0L]
  signature_matrix <- do.call(rbind, lapply(marker_sets, function(markers) {
    Matrix::colMeans(expression[markers, , drop = FALSE])
  }))
  rownames(signature_matrix) <- names(marker_sets)

  bind_rows(lapply(names(marker_sets), function(annotation) {
    assigned <- labels == annotation
    assigned[is.na(assigned)] <- FALSE
    comparator <- !assigned & !is.na(labels)
    values <- signature_matrix[annotation, ]
    all_sd <- sd(values[assigned | comparator])
    tibble(
      compartment = compartment,
      strategy = strategy,
      strategy_label = unname(strategy_labels[strategy]),
      annotation = annotation,
      markers = paste(marker_sets[[annotation]], collapse = "; "),
      n_assigned = sum(assigned),
      n_comparator = sum(comparator),
      mean_matching_signature_assigned = if (any(assigned)) {
        mean(values[assigned])
      } else {
        NA_real_
      },
      mean_matching_signature_other = if (any(comparator)) {
        mean(values[comparator])
      } else {
        NA_real_
      },
      delta_matching_signature = mean_matching_signature_assigned -
        mean_matching_signature_other,
      standardized_delta = if (is.finite(all_sd) && all_sd > 0) {
        delta_matching_signature / all_sd
      } else {
        NA_real_
      }
    )
  }))
}

plot_marker_dot <- function(dot_data, title) {
  ggplot(
    dot_data,
    aes(x = gene, y = annotation, size = percent_expressed, color = scaled_mean_expression)
  ) +
    geom_point() +
    facet_wrap(~strategy_label, ncol = 1) +
    scale_size(
      name = "Percent expressed",
      range = c(1.5, 7),
      limits = c(0, 100),
      breaks = c(25, 50, 75, 100)
    ) +
    scale_color_gradient2(
      name = "Scaled mean\nexpression",
      low = "#3B6FB6",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-2, 2),
      oob = scales::squish
    ) +
    labs(title = title, x = NULL, y = NULL) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      strip.background = element_blank(),
      strip.text = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 10),
      legend.position = "right"
    )
}

message("Reading candidate seven-donor SFE")
sfe7 <- readRDS(annotated_path)
if (!file.exists(joint_candidate_path)) {
  saveRDS(sfe7, joint_candidate_path)
}

genes <- rownames(sfe7)
availability <- bind_rows(
  marker_availability(major_markers, genes, "Major populations"),
  marker_availability(fib_markers, genes, "Fibroblast subtypes"),
  marker_availability(ec_markers, genes, "Endothelial subtypes")
)
write.csv(
  availability,
  file.path(output_dir, "annotation_marker_availability.csv"),
  row.names = FALSE
)
stopifnot(all(availability$available))

query <- as.character(sfe7$SampleId) == "ExcludedValidation"
expression <- assay(sfe7, "logcounts")[, query, drop = FALSE]
counts <- assay(sfe7, "counts")[, query, drop = FALSE]

major_dot <- bind_rows(
  make_dot_summary(
    expression,
    counts,
    as.character(sfe7$lv1_anno_joint[query]),
    major_markers,
    "joint_harmony",
    major_levels
  ),
  make_dot_summary(
    expression,
    counts,
    as.character(sfe7$lv1_anno_mapped[query]),
    major_markers,
    "reference_mapping",
    major_levels
  )
)
fib_dot <- bind_rows(
  make_dot_summary(
    expression,
    counts,
    as.character(sfe7$lv2_anno_joint[query]),
    fib_markers,
    "joint_harmony",
    fib_levels
  ),
  make_dot_summary(
    expression,
    counts,
    as.character(sfe7$lv2_anno_mapped[query]),
    fib_markers,
    "reference_mapping",
    fib_levels
  )
)
ec_dot <- bind_rows(
  make_dot_summary(
    expression,
    counts,
    as.character(sfe7$lv2_anno_joint[query]),
    ec_markers,
    "joint_harmony",
    ec_levels
  ),
  make_dot_summary(
    expression,
    counts,
    as.character(sfe7$lv2_anno_mapped[query]),
    ec_markers,
    "reference_mapping",
    ec_levels
  )
)

write.csv(
  major_dot,
  file.path(output_dir, "ssc230_major_marker_dotplot_data.csv"),
  row.names = FALSE
)
write.csv(
  fib_dot,
  file.path(output_dir, "ssc230_fibroblast_marker_dotplot_data.csv"),
  row.names = FALSE
)
write.csv(
  ec_dot,
  file.path(output_dir, "ssc230_endothelial_marker_dotplot_data.csv"),
  row.names = FALSE
)

save_plot_set(
  plot_marker_dot(major_dot, "ExcludedValidation major-population marker QC"),
  "ssc230_major_marker_qc",
  width = 11,
  height = 8
)
save_plot_set(
  plot_marker_dot(fib_dot, "ExcludedValidation fibroblast-subtype marker QC"),
  "ssc230_fibroblast_marker_qc",
  width = 10,
  height = 7.5
)
save_plot_set(
  plot_marker_dot(ec_dot, "ExcludedValidation endothelial-subtype marker QC"),
  "ssc230_endothelial_marker_qc",
  width = 9,
  height = 6.5
)

specificity <- bind_rows(
  signature_specificity(
    expression,
    as.character(sfe7$lv1_anno_joint[query]),
    major_markers,
    "joint_harmony",
    "Major populations"
  ),
  signature_specificity(
    expression,
    as.character(sfe7$lv1_anno_mapped[query]),
    major_markers,
    "reference_mapping",
    "Major populations"
  ),
  signature_specificity(
    expression,
    as.character(sfe7$lv2_anno_joint[query]),
    fib_markers,
    "joint_harmony",
    "Fibroblast subtypes"
  ),
  signature_specificity(
    expression,
    as.character(sfe7$lv2_anno_mapped[query]),
    fib_markers,
    "reference_mapping",
    "Fibroblast subtypes"
  ),
  signature_specificity(
    expression,
    as.character(sfe7$lv2_anno_joint[query]),
    ec_markers,
    "joint_harmony",
    "Endothelial subtypes"
  ),
  signature_specificity(
    expression,
    as.character(sfe7$lv2_anno_mapped[query]),
    ec_markers,
    "reference_mapping",
    "Endothelial subtypes"
  )
)

specificity_summary <- specificity %>%
  filter(n_assigned > 0, is.finite(standardized_delta)) %>%
  group_by(compartment, strategy, strategy_label) %>%
  summarise(
    n_labels_represented = n(),
    n_cells_assigned = sum(n_assigned),
    weighted_standardized_delta = weighted.mean(
      standardized_delta,
      w = n_assigned
    ),
    proportion_labels_positive = mean(delta_matching_signature > 0),
    .groups = "drop"
  )

write.csv(
  specificity,
  file.path(output_dir, "ssc230_annotation_marker_specificity.csv"),
  row.names = FALSE
)
write.csv(
  specificity_summary,
  file.path(output_dir, "ssc230_annotation_marker_specificity_summary.csv"),
  row.names = FALSE
)

comparison <- read.csv(
  file.path(output_dir, "annotation_strategy_comparison_summary.csv"),
  check.names = FALSE
)
metric_value <- setNames(comparison$value, comparison$metric)
major_comparable <- unname(metric_value[
  "ExcludedValidation major-label agreement: joint versus mapping"
]) >= 0.85
subtype_refit_stable <- unname(metric_value[
  "Six-donor Fib/EC subtype retention in joint clustering"
]) >= 0.75

# The original six-donor fine annotations were manually curated after ResolVI.
# Harmony is an explicit sensitivity analysis because that trained ResolVI
# model/H5AD is unavailable. Preserve the curated labels and map only ExcludedValidation
# when broad compartments are concordant and a Harmony refit is not subtype-stable.
recommended_strategy <- if (major_comparable && !subtype_refit_stable) {
  "reference_mapping_with_joint_harmony_sensitivity"
} else if (major_comparable) {
  "reference_mapping"
} else {
  "manual_review_required"
}

if (recommended_strategy == "manual_review_required") {
  stop(
    "Major-population agreement was below the prespecified 0.85 threshold; ",
    "final labels were not overwritten."
  )
}

decision <- tibble(
  recommended_strategy = recommended_strategy,
  ssc230_major_agreement_at_least_0_85 = major_comparable,
  six_donor_subtype_retention_at_least_0_75 = subtype_refit_stable,
  rationale = paste(
    "Reference mapping is used for manuscript-facing labels because broad",
    "ExcludedValidation compartments are concordant between strategies, while replacing",
    "the original ResolVI/manual fine annotations with Harmony would alter",
    "more than half of the six-donor Fib/EC subtype calls. Joint Harmony is",
    "retained as a seven-donor sensitivity analysis and embedding."
  )
)
write.csv(
  decision,
  file.path(output_dir, "annotation_strategy_decision.csv"),
  row.names = FALSE
)

message("Finalizing the SFE with reference-mapped manuscript labels")
sfe7$lv1_anno_final <- factor(
  as.character(sfe7$lv1_anno_mapped),
  levels = major_levels
)
sfe7$lv2_anno_final <- factor(
  as.character(sfe7$lv2_anno_mapped),
  levels = c(fib_levels, ec_levels)
)
sfe7$lv1_anno_new_text <- sfe7$lv1_anno_final
sfe7$lv2_anno_new <- sfe7$lv2_anno_final
sfe7$lv3_anno <- sfe7$lv2_anno_final
sfe7$annotation_strategy_final <- recommended_strategy

metadata(sfe7)$annotation_comparison$marker_qc <- list(
  marker_sets = list(
    major = major_markers,
    fibroblast = fib_markers,
    endothelial = ec_markers
  ),
  specificity_summary = specificity_summary,
  decision = decision
)
metadata(sfe7)$annotation_comparison$selected_strategy <- recommended_strategy
metadata(sfe7)$annotation_comparison$decision <- decision

saveRDS(sfe7, annotated_path)

final_counts <- as.data.frame(table(
  donor = as.character(sfe7$SampleId),
  lv1_anno_final = as.character(sfe7$lv1_anno_final)
)) %>%
  filter(Freq > 0) %>%
  rename(n_cells = Freq)
write.csv(
  final_counts,
  file.path(output_dir, "seven_donor_final_major_label_counts.csv"),
  row.names = FALSE
)

cat("\nMarker specificity summary\n")
print(specificity_summary)
cat("\nFinal annotation decision\n")
print(decision)
cat("\nSaved final annotated SFE:", normalizePath(annotated_path, winslash = "/"), "\n")
