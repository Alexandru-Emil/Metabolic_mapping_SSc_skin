#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SpatialFeatureExperiment)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(svglite)
  library(jsonlite)
})

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "Usage: 03_attach_joint_resolvi_and_marker_qc.R ",
    "<input_sfe.rds> <annotation_dir> <output_dir>"
  )
}

input_sfe_path <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
annotation_dir <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[3]], winslash = "/", mustWork = FALSE)
plot_dir <- file.path(output_dir, "annotation_qc_plots")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

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
  ACTA2_EC = "ACTA2",
  HEY1_EC = c("HEY1", "GJA5"),
  lymphatic_EC = c("LYVE1", "PDPN")
)
strategy_labels <- c(
  joint_resolvi = "Joint-ResolVI latent mapping",
  prior_mapping = "Prior six-donor reference mapping"
)

save_plot_set <- function(plot, stem, width, height) {
  ggsave(
    file.path(plot_dir, paste0(stem, ".png")),
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  ggsave(
    file.path(plot_dir, paste0(stem, ".pdf")),
    plot,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf,
    bg = "white"
  )
  ggsave(
    file.path(plot_dir, paste0(stem, ".svg")),
    plot,
    width = width,
    height = height,
    units = "in",
    device = svglite::svglite,
    bg = "white"
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

  bind_rows(lapply(levels(droplevels(labels)), function(label) {
    index <- labels == label
    tibble(
      annotation = label,
      gene = marker_order,
      mean_expression = Matrix::rowMeans(
        expression[, index, drop = FALSE]
      ),
      percent_expressed = 100 * Matrix::rowMeans(
        counts[, index, drop = FALSE] > 0
      ),
      n_cells = sum(index)
    )
  })) %>%
    group_by(gene) %>%
    mutate(
      scaled_mean_expression = if (
        n() > 1L && is.finite(sd(mean_expression)) &&
          sd(mean_expression) > 0
      ) {
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
  marker_sets <- lapply(marker_sets, function(markers) {
    markers[markers %in% rownames(expression)]
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
    combined_sd <- sd(values[assigned | comparator])
    assigned_mean <- if (any(assigned)) mean(values[assigned]) else NA_real_
    comparator_mean <- if (any(comparator)) {
      mean(values[comparator])
    } else {
      NA_real_
    }
    delta <- assigned_mean - comparator_mean
    tibble(
      compartment = compartment,
      strategy = strategy,
      strategy_label = unname(strategy_labels[strategy]),
      annotation = annotation,
      markers = paste(marker_sets[[annotation]], collapse = "; "),
      n_assigned = sum(assigned),
      n_comparator = sum(comparator),
      mean_matching_signature_assigned = assigned_mean,
      mean_matching_signature_other = comparator_mean,
      delta_matching_signature = delta,
      standardized_delta = if (
        is.finite(combined_sd) && combined_sd > 0
      ) {
        delta / combined_sd
      } else {
        NA_real_
      }
    )
  }))
}

plot_marker_dot <- function(data, title) {
  ggplot(
    data,
    aes(
      x = gene,
      y = annotation,
      size = percent_expressed,
      color = scaled_mean_expression
    )
  ) +
    geom_point() +
    facet_wrap(~strategy_label, ncol = 1, scales = "free_y") +
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
      strip.text = element_text(size = 11, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 10),
      legend.position = "right"
    )
}

message("Reading local seven-donor SFE")
sfe <- readRDS(input_sfe_path)
stopifnot(
  inherits(sfe, "SpatialFeatureExperiment"),
  nrow(sfe) == 5099L,
  length(unique(as.character(sfe$SampleId))) == 7L,
  !anyDuplicated(colnames(sfe))
)

message("Reading joint ResolVI annotations and latent")
annotations <- read.csv(
  gzfile(file.path(annotation_dir, "joint_resolvi_cell_annotations.csv.gz")),
  check.names = FALSE,
  na.strings = c("", "NA", "<NA>")
)
latent <- read.csv(
  gzfile(file.path(annotation_dir, "joint_resolvi_latent.csv.gz")),
  check.names = FALSE
)
annotation_index <- match(colnames(sfe), annotations$cell)
latent_index <- match(colnames(sfe), latent$cell)
stopifnot(
  !anyNA(annotation_index),
  !anyNA(latent_index),
  !anyDuplicated(annotations$cell),
  !anyDuplicated(latent$cell)
)
annotations <- annotations[annotation_index, , drop = FALSE]
latent <- latent[latent_index, , drop = FALSE]
stopifnot(
  identical(annotations$cell, colnames(sfe)),
  identical(latent$cell, colnames(sfe)),
  identical(as.character(annotations$SampleId), as.character(sfe$SampleId))
)

input_reference_lv1 <- as.character(sfe$lv1_anno_final)
input_reference_lv2 <- as.character(sfe$lv2_anno_final)
sfe$lv1_anno_reference_mapped <- factor(
  input_reference_lv1,
  levels = major_levels
)
sfe$lv2_anno_reference_mapped <- factor(
  input_reference_lv2,
  levels = c(fib_levels, ec_levels)
)
if ("lv1_anno_joint" %in% colnames(colData(sfe))) {
  sfe$lv1_anno_harmony_sensitivity <- sfe$lv1_anno_joint
}
if ("lv2_anno_joint" %in% colnames(colData(sfe))) {
  sfe$lv2_anno_harmony_sensitivity <- sfe$lv2_anno_joint
}

sfe$lv1_anno_resolvi_mapped <- factor(
  annotations$joint_resolvi_lv1,
  levels = major_levels
)
sfe$lv2_anno_resolvi_mapped <- factor(
  annotations$joint_resolvi_lv2,
  levels = c(fib_levels, ec_levels)
)
sfe$lv1_anno_joint_resolvi <- sfe$lv1_anno_resolvi_mapped
sfe$lv2_anno_joint_resolvi <- sfe$lv2_anno_resolvi_mapped
sfe$lv1_anno_cluster_majority_resolvi <- factor(
  annotations$cluster_majority_lv1_sensitivity,
  levels = major_levels
)
sfe$lv2_anno_cluster_majority_resolvi <- factor(
  annotations$cluster_majority_lv2_sensitivity,
  levels = c(fib_levels, ec_levels)
)
sfe$resolvi_lv1_mapping_confidence <- as.numeric(
  annotations$joint_resolvi_lv1_mapping_confidence
)
sfe$resolvi_lv2_mapping_confidence <- as.numeric(
  annotations$joint_resolvi_lv2_mapping_confidence
)
sfe$resolvi_major_cluster <- factor(
  annotations$joint_resolvi_major_cluster
)
sfe$lv1_anno_final <- sfe$lv1_anno_resolvi_mapped
sfe$lv2_anno_final <- sfe$lv2_anno_resolvi_mapped
sfe$lv1_anno_new_text <- sfe$lv1_anno_final
sfe$lv2_anno_new <- sfe$lv2_anno_final
sfe$lv3_anno <- sfe$lv2_anno_final
sfe$annotation_strategy_final <- paste0(
  "curated_six_donor_plus_ssc230_joint_resolvi_latent_mapping"
)

reference_donor <- as.character(sfe$SampleId) != "ExcludedValidation"
stopifnot(
  identical(
    as.character(sfe$lv1_anno_final[reference_donor]),
    input_reference_lv1[reference_donor]
  ),
  identical(
    as.character(sfe$lv2_anno_final[reference_donor]),
    input_reference_lv2[reference_donor]
  ),
  !anyNA(sfe$lv1_anno_final)
)

latent_matrix <- as.matrix(latent[, setdiff(colnames(latent), "cell")])
rownames(latent_matrix) <- latent$cell
storage.mode(latent_matrix) <- "double"
umap_matrix <- as.matrix(annotations[, c(
  "joint_resolvi_umap1",
  "joint_resolvi_umap2"
)])
rownames(umap_matrix) <- annotations$cell
colnames(umap_matrix) <- c("UMAP1", "UMAP2")
storage.mode(umap_matrix) <- "double"
stopifnot(
  all(is.finite(latent_matrix)),
  all(is.finite(umap_matrix)),
  ncol(latent_matrix) == 30L
)
reducedDim(sfe, "X_resolvi") <- latent_matrix
reducedDim(sfe, "umap_joint_resolvi") <- umap_matrix

annotation_manifest <- jsonlite::read_json(
  file.path(annotation_dir, "joint_resolvi_annotation_manifest.json"),
  simplifyVector = TRUE
)
metadata(sfe)$joint_seven_donor_resolvi <- list(
  model_annotation_manifest = annotation_manifest,
  final_annotation = paste(
    "The six curated donors retain their input labels exactly. ExcludedValidation labels",
    "are predicted from the joint ResolVI latent with a k-nearest-neighbor",
    "mapper tuned by leave-one-donor-out label-transfer validation among the",
    "six reference donors. All seven donors contributed to unsupervised",
    "ResolVI training and latent geometry."
  ),
  cluster_sensitivity = paste(
    "Leiden cluster-majority labels are retained only in the",
    "lv1_anno_cluster_majority_resolvi and",
    "lv2_anno_cluster_majority_resolvi columns. They are not final labels",
    "because mixed clusters collapse rare major lineages and fine subtypes."
  ),
  final_label_columns = c("lv1_anno_final", "lv2_anno_final"),
  preserved_prior_mapping_columns = c(
    "lv1_anno_reference_mapped",
    "lv2_anno_reference_mapped"
  ),
  final_strategy = unique(as.character(sfe$annotation_strategy_final))
)

message("Saving joint ResolVI SFE")
final_path <- file.path(output_dir, "sfe_seven_donors_annotated_final.rds")
explicit_path <- file.path(
  output_dir,
  "sfe_seven_donors_annotated_joint_resolvi.rds"
)
saveRDS(sfe, final_path)
saveRDS(sfe, explicit_path)

join_audit <- as.data.frame(colData(sfe), check.names = FALSE) %>%
  transmute(
    donor = as.character(SampleId),
    valid_key = !is.na(KEY) & nzchar(as.character(KEY)) &
      as.character(KEY) != "not found_not found",
    matched_imc = as.logical(imc_matched)
  ) %>%
  group_by(donor) %>%
  summarise(
    n_xenium_cells = n(),
    n_valid_keys = sum(valid_key, na.rm = TRUE),
    n_matched_imc = sum(matched_imc, na.rm = TRUE),
    matched_imc_percent = round(100 * n_matched_imc / n_xenium_cells, 1),
    .groups = "drop"
  )
write.csv(
  join_audit,
  file.path(output_dir, "seven_donor_join_audit.csv"),
  row.names = FALSE
)

assay_audit <- tibble(
  assay = assayNames(sfe),
  matrix_class = vapply(
    assayNames(sfe),
    function(name) paste(class(assay(sfe, name)), collapse = "/"),
    character(1)
  ),
  n_features = nrow(sfe),
  n_cells = ncol(sfe),
  nonzero = vapply(
    assayNames(sfe),
    function(name) Matrix::nnzero(assay(sfe, name)),
    numeric(1)
  )
)
write.csv(
  assay_audit,
  file.path(output_dir, "seven_donor_assay_audit.csv"),
  row.names = FALSE
)

major_counts <- as.data.frame(table(
  donor = as.character(sfe$SampleId),
  lv1_anno_final = as.character(sfe$lv1_anno_final)
)) %>%
  filter(Freq > 0) %>%
  rename(n_cells = Freq)
write.csv(
  major_counts,
  file.path(output_dir, "seven_donor_final_major_label_counts.csv"),
  row.names = FALSE
)

query <- as.character(sfe$SampleId) == "ExcludedValidation"
reference <- !query
comparison <- read.csv(
  file.path(annotation_dir, "joint_resolvi_annotation_summary.csv"),
  check.names = FALSE
)
write.csv(
  comparison,
  file.path(output_dir, "joint_resolvi_annotation_summary.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(
    recommended_strategy = paste0(
      "curated_six_donor_plus_ssc230_joint_resolvi_latent_mapping"
    ),
    rationale = paste(
      "All seven donors contribute to the spatial generative model and latent",
      "geometry. The six curated donors remain unchanged; only ExcludedValidation is",
      "mapped in the joint latent with a sample-aware classifier. Joint",
      "cluster-majority labels are sensitivity outputs because they merge",
      "rare lineages and do not retain the curated fine-subtype vocabulary."
    )
  ),
  file.path(output_dir, "annotation_strategy_decision.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(table(
    joint_resolvi = as.character(sfe$lv1_anno_final[query]),
    prior_mapping = as.character(sfe$lv1_anno_reference_mapped[query])
  )),
  file.path(output_dir, "ssc230_major_joint_resolvi_vs_mapping.csv"),
  row.names = FALSE
)
write.csv(
  as.data.frame(table(
    joint_resolvi = as.character(sfe$lv2_anno_final[query]),
    prior_mapping = as.character(sfe$lv2_anno_reference_mapped[query])
  )),
  file.path(output_dir, "ssc230_subtype_joint_resolvi_vs_mapping.csv"),
  row.names = FALSE
)
write.csv(
  as.data.frame(table(
    cluster_majority_resolvi = as.character(
      sfe$lv1_anno_cluster_majority_resolvi[reference]
    ),
    curated_reference = as.character(
      sfe$lv1_anno_reference_mapped[reference]
    )
  )),
  file.path(output_dir, "six_donor_major_label_retention_crosstab.csv"),
  row.names = FALSE
)
write.csv(
  as.data.frame(table(
    cluster_majority_resolvi = as.character(
      sfe$lv2_anno_cluster_majority_resolvi[reference]
    ),
    curated_reference = as.character(
      sfe$lv2_anno_reference_mapped[reference]
    )
  )),
  file.path(output_dir, "six_donor_subtype_label_retention_crosstab.csv"),
  row.names = FALSE
)

availability <- bind_rows(
  tibble(
    compartment = "Major populations",
    annotation = rep(names(major_markers), lengths(major_markers)),
    gene = unlist(major_markers, use.names = FALSE)
  ),
  tibble(
    compartment = "Fibroblast subtypes",
    annotation = rep(names(fib_markers), lengths(fib_markers)),
    gene = unlist(fib_markers, use.names = FALSE)
  ),
  tibble(
    compartment = "Endothelial subtypes",
    annotation = rep(names(ec_markers), lengths(ec_markers)),
    gene = unlist(ec_markers, use.names = FALSE)
  )
) %>%
  mutate(available = gene %in% rownames(sfe))
write.csv(
  availability,
  file.path(output_dir, "annotation_marker_availability.csv"),
  row.names = FALSE
)
stopifnot(all(availability$available))

message("Generating ExcludedValidation marker QC")
expression <- assay(sfe, "logcounts")[, query, drop = FALSE]
counts <- assay(sfe, "counts")[, query, drop = FALSE]

major_dot <- bind_rows(
  make_dot_summary(
    expression,
    counts,
    as.character(sfe$lv1_anno_final[query]),
    major_markers,
    "joint_resolvi",
    major_levels
  ),
  make_dot_summary(
    expression,
    counts,
    as.character(sfe$lv1_anno_reference_mapped[query]),
    major_markers,
    "prior_mapping",
    major_levels
  )
)
fib_dot <- bind_rows(
  make_dot_summary(
    expression,
    counts,
    as.character(sfe$lv2_anno_final[query]),
    fib_markers,
    "joint_resolvi",
    fib_levels
  ),
  make_dot_summary(
    expression,
    counts,
    as.character(sfe$lv2_anno_reference_mapped[query]),
    fib_markers,
    "prior_mapping",
    fib_levels
  )
)
ec_dot <- bind_rows(
  make_dot_summary(
    expression,
    counts,
    as.character(sfe$lv2_anno_final[query]),
    ec_markers,
    "joint_resolvi",
    ec_levels
  ),
  make_dot_summary(
    expression,
    counts,
    as.character(sfe$lv2_anno_reference_mapped[query]),
    ec_markers,
    "prior_mapping",
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
  11,
  8
)
save_plot_set(
  plot_marker_dot(fib_dot, "ExcludedValidation fibroblast-subtype marker QC"),
  "ssc230_fibroblast_marker_qc",
  10,
  7.5
)
save_plot_set(
  plot_marker_dot(ec_dot, "ExcludedValidation endothelial-subtype marker QC"),
  "ssc230_endothelial_marker_qc",
  9,
  6.5
)

specificity <- bind_rows(
  signature_specificity(
    expression,
    as.character(sfe$lv1_anno_final[query]),
    major_markers,
    "joint_resolvi",
    "Major populations"
  ),
  signature_specificity(
    expression,
    as.character(sfe$lv1_anno_reference_mapped[query]),
    major_markers,
    "prior_mapping",
    "Major populations"
  ),
  signature_specificity(
    expression,
    as.character(sfe$lv2_anno_final[query]),
    fib_markers,
    "joint_resolvi",
    "Fibroblast subtypes"
  ),
  signature_specificity(
    expression,
    as.character(sfe$lv2_anno_reference_mapped[query]),
    fib_markers,
    "prior_mapping",
    "Fibroblast subtypes"
  ),
  signature_specificity(
    expression,
    as.character(sfe$lv2_anno_final[query]),
    ec_markers,
    "joint_resolvi",
    "Endothelial subtypes"
  ),
  signature_specificity(
    expression,
    as.character(sfe$lv2_anno_reference_mapped[query]),
    ec_markers,
    "prior_mapping",
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

cat("\nJoint ResolVI annotation summary\n")
print(comparison)
cat("\nFinal major-label counts\n")
print(major_counts)
cat("\nSaved:", normalizePath(final_path, winslash = "/"), "\n")
