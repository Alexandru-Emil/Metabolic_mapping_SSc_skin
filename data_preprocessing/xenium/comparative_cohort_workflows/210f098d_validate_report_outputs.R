#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SpatialFeatureExperiment)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(jsonlite)
  library(dplyr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) {
  stop("Usage: validate_report_outputs.R [report_dir]")
}

script_argument <- grep("^--file=", commandArgs(), value = TRUE)
stopifnot(length(script_argument) == 1L)
report_dir <- normalizePath(
  if (length(args) == 1L) {
    args[[1]]
  } else {
    dirname(dirname(sub("^--file=", "", script_argument)))
  },
  winslash = "/",
  mustWork = TRUE
)
output_dir <- file.path(report_dir, "outputs")
panel_dir <- file.path(output_dir, "seven_donor_manuscript_panels")
model_dir <- file.path(report_dir, "model_artifacts", "resolvi_run")
input_sfe_path <- file.path(
  report_dir,
  "inputs",
  "sfe_seven_donors_annotated_reference_mapped_20260729.rds"
)

checks <- list()
add_check <- function(name, passed, detail) {
  checks[[length(checks) + 1L]] <<- tibble(
    check = name,
    passed = isTRUE(passed),
    detail = as.character(detail)
  )
}

required_files <- c(
  file.path(
    report_dir,
    "Seven_donor_IMC_Xenium_joint_ResolVI_reanalysis.Rmd"
  ),
  file.path(
    report_dir,
    "Seven_donor_IMC_Xenium_joint_ResolVI_reanalysis.html"
  ),
  input_sfe_path,
  file.path(output_dir, "sfe_seven_donors_annotated_final.rds"),
  file.path(output_dir, "joint_resolvi_cell_annotations.csv.gz"),
  file.path(output_dir, "joint_resolvi_latent.csv.gz"),
  file.path(output_dir, "resolvi_training_manifest.json"),
  file.path(output_dir, "joint_resolvi_annotation_manifest.json"),
  file.path(output_dir, "annotation_strategy_decision.csv"),
  file.path(output_dir, "major_knn_lodo_grid.csv"),
  file.path(output_dir, "fibroblast_subtype_knn_lodo_grid.csv"),
  file.path(output_dir, "endothelial_subtype_knn_lodo_grid.csv"),
  file.path(
    output_dir,
    "seven_donor_scores",
    "seven_donor_score_audit.csv"
  ),
  file.path(
    panel_dir,
    "seven_donor_figure_export_manifest.csv"
  ),
  file.path(
    model_dir,
    "validation",
    "resolvi_validation.json"
  ),
  file.path(model_dir, "seven_donor_resolvi.h5ad"),
  file.path(model_dir, "model", "model.pt")
)
missing_files <- required_files[!file.exists(required_files)]
add_check(
  "Required report/model files",
  !length(missing_files),
  if (length(missing_files)) {
    paste(missing_files, collapse = "; ")
  } else {
    paste(length(required_files), "files present")
  }
)

html_path <- file.path(
  report_dir,
  "Seven_donor_IMC_Xenium_joint_ResolVI_reanalysis.html"
)
html_size <- if (file.exists(html_path)) file.info(html_path)$size else 0
add_check(
  "Rendered HTML is non-empty",
  html_size > 100000L,
  paste(html_size, "bytes")
)

sfe_path <- file.path(output_dir, "sfe_seven_donors_annotated_final.rds")
if (file.exists(sfe_path)) {
  sfe <- readRDS(sfe_path)
  expected_donors <- c(
    "Validation1" = 2872L,
    "Validation2" = 3934L,
    "ExcludedValidation" = 1880L,
    "Validation3" = 1775L,
    "Validation4" = 3413L,
    "Validation5" = 4174L,
    "Validation6" = 3223L
  )
  observed_donors <- table(as.character(sfe$SampleId))
  add_check(
    "Final SFE class and dimensions",
    inherits(sfe, "SpatialFeatureExperiment") &&
      identical(dim(sfe), c(5099L, 21271L)),
    paste(class(sfe)[1], paste(dim(sfe), collapse = " x "))
  )
  add_check(
    "Final SFE donor counts",
    identical(
      as.integer(observed_donors[names(expected_donors)]),
      unname(expected_donors)
    ),
    paste(
      names(expected_donors),
      as.integer(observed_donors[names(expected_donors)]),
      collapse = "; "
    )
  )
  add_check(
    "ResolVI reduced dimensions",
    all(c("X_resolvi", "umap_joint_resolvi") %in% reducedDimNames(sfe)) &&
      ncol(reducedDim(sfe, "X_resolvi")) == 30L &&
      all(is.finite(reducedDim(sfe, "X_resolvi"))) &&
      all(is.finite(reducedDim(sfe, "umap_joint_resolvi"))),
    paste(reducedDimNames(sfe), collapse = ", ")
  )
  required_metadata <- c(
    "lv1_anno_final",
    "lv2_anno_final",
    "lv1_anno_reference_mapped",
    "lv2_anno_reference_mapped",
    "lv1_anno_resolvi_mapped",
    "lv2_anno_resolvi_mapped",
    "lv1_anno_cluster_majority_resolvi",
    "lv2_anno_cluster_majority_resolvi",
    "resolvi_lv1_mapping_confidence",
    "resolvi_lv2_mapping_confidence",
    "annotation_strategy_final"
  )
  expected_strategy <- paste0(
    "curated_six_donor_plus_ssc230_joint_resolvi_latent_mapping"
  )
  add_check(
    "Final, preserved, and sensitivity annotation columns",
    all(required_metadata %in% colnames(colData(sfe))) &&
      !anyNA(sfe$lv1_anno_final) &&
      all(sfe$annotation_strategy_final == expected_strategy),
    paste(intersect(required_metadata, colnames(colData(sfe))), collapse = ", ")
  )
  if (file.exists(input_sfe_path)) {
    input_sfe <- readRDS(input_sfe_path)
    reference_cells <- as.character(input_sfe$SampleId) != "ExcludedValidation"
    query_cells <- !reference_cells
    add_check(
      "Input cell/gene order preserved",
      identical(rownames(input_sfe), rownames(sfe)) &&
        identical(colnames(input_sfe), colnames(sfe)),
      paste(nrow(sfe), "genes and", ncol(sfe), "cells")
    )
    add_check(
      "Raw counts preserved exactly",
      identical(assay(input_sfe, "counts"), assay(sfe, "counts")),
      paste(
        "sum =",
        format(sum(assay(sfe, "counts")), scientific = FALSE),
        "; nonzero =",
        Matrix::nnzero(assay(sfe, "counts"))
      )
    )
    add_check(
      "Spatial coordinates preserved exactly",
      identical(spatialCoords(input_sfe), spatialCoords(sfe)),
      paste(dim(spatialCoords(sfe)), collapse = " x ")
    )
    geometry_names <- colGeometryNames(input_sfe)
    geometries_equal <- identical(geometry_names, colGeometryNames(sfe)) &&
      all(vapply(
        geometry_names,
        function(name) {
          isTRUE(all.equal(
            colGeometry(input_sfe, name),
            colGeometry(sfe, name),
            check.attributes = TRUE
          ))
        },
        logical(1)
      ))
    add_check(
      "Cell centroids and boundaries preserved",
      geometries_equal,
      paste(geometry_names, collapse = ", ")
    )
    add_check(
      "Prior reference-mapped labels preserved",
      identical(
        as.character(sfe$lv1_anno_reference_mapped),
        as.character(input_sfe$lv1_anno_final)
      ) &&
        identical(
          as.character(sfe$lv2_anno_reference_mapped),
          as.character(input_sfe$lv2_anno_final)
        ),
      "Both input label columns retained for all 21,271 cells"
    )
    add_check(
      "Six curated donors are not relabelled",
      identical(
        as.character(sfe$lv1_anno_final[reference_cells]),
        as.character(input_sfe$lv1_anno_final[reference_cells])
      ) &&
        identical(
          as.character(sfe$lv2_anno_final[reference_cells]),
          as.character(input_sfe$lv2_anno_final[reference_cells])
        ),
      paste(sum(reference_cells), "reference cells compared")
    )
    query_lineage <- query_cells &
      as.character(sfe$lv1_anno_final) %in% c("fibroblast", "endothelial")
    add_check(
      "ExcludedValidation ResolVI mapping is complete",
      all(is.finite(sfe$resolvi_lv1_mapping_confidence[query_cells])) &&
        all(is.finite(
          sfe$resolvi_lv2_mapping_confidence[query_lineage]
        )),
      paste(
        sum(query_cells),
        "major and",
        sum(query_lineage),
        "Fib/EC subtype predictions"
      )
    )
    major_levels <- c(
      "epithelial", "fibroblast", "pericyte", "endothelial", "muscle",
      "telocyte", "Schwann cell", "myeloid cell", "lymphocyte"
    )
    subtype_levels <- c(
      "papillary_Fib", "COMP_Fib", "COL8A1_Fib", "PI16_Fib",
      "CXCL12_Fib", "CCL19_Fib", "NGFR_Fib", "COCH_Fib",
      "cycling_EC", "ACKR1_EC", "EPC", "ACTA2_EC", "HEY1_EC",
      "lymphatic_EC"
    )
    add_check(
      "Final annotation retains complete curated vocabulary",
      all(major_levels %in% unique(as.character(sfe$lv1_anno_final))) &&
        all(subtype_levels %in% unique(as.character(sfe$lv2_anno_final))),
      paste(length(major_levels), "major and", length(subtype_levels), "subtypes")
    )
    add_check(
      "Cluster-majority labels remain sensitivity-only",
      any(
        as.character(sfe$lv1_anno_cluster_majority_resolvi) !=
          as.character(sfe$lv1_anno_final),
        na.rm = TRUE
      ) &&
        identical(
          as.character(sfe$lv1_anno_resolvi_mapped),
          as.character(sfe$lv1_anno_final)
        ),
      "Final labels equal latent mapping, not cluster-majority calls"
    )
    rm(input_sfe)
  }
}

training_manifest_path <- file.path(
  output_dir,
  "resolvi_training_manifest.json"
)
annotation_manifest_path <- file.path(
  output_dir,
  "joint_resolvi_annotation_manifest.json"
)
validation_path <- file.path(
  model_dir,
  "validation",
  "resolvi_validation.json"
)
if (
  file.exists(training_manifest_path) &&
    file.exists(annotation_manifest_path) &&
    file.exists(validation_path)
) {
  training_manifest <- read_json(training_manifest_path, simplifyVector = TRUE)
  annotation_manifest <- read_json(
    annotation_manifest_path,
    simplifyVector = TRUE
  )
  validation <- read_json(validation_path, simplifyVector = TRUE)
  add_check(
    "Training and model validation status",
    identical(training_manifest$status, "complete") &&
      identical(validation$status, "passed"),
    paste(training_manifest$status, validation$status, sep = " / ")
  )
  add_check(
    "Annotation used the validated model H5AD",
    identical(
      annotation_manifest$model_h5ad_sha256,
      validation$output_h5ad_sha256
    ),
    annotation_manifest$model_h5ad_sha256
  )
  add_check(
    "Annotation manifest preserves six-donor reference",
    isTRUE(annotation_manifest$six_donor_reference_labels_preserved) &&
      identical(
        annotation_manifest$joint_leiden_role,
        "sensitivity analysis and visualization only"
      ),
    annotation_manifest$final_annotation_strategy
  )
}

strategy_path <- file.path(output_dir, "annotation_strategy_decision.csv")
if (file.exists(strategy_path)) {
  strategy <- read.csv(strategy_path, check.names = FALSE)
  add_check(
    "Recorded final annotation strategy",
    nrow(strategy) == 1L &&
      identical(
        strategy$recommended_strategy[[1]],
        "curated_six_donor_plus_ssc230_joint_resolvi_latent_mapping"
      ),
    strategy$recommended_strategy[[1]]
  )
}

panel_manifest_path <- file.path(
  panel_dir,
  "seven_donor_figure_export_manifest.csv"
)
if (file.exists(panel_manifest_path)) {
  panel_manifest <- read.csv(panel_manifest_path, check.names = FALSE)
  expected_panels <- unlist(lapply(
    panel_manifest$file_stem,
    function(stem) file.path(panel_dir, paste0(stem, c(".png", ".pdf", ".svg")))
  ))
  panel_sizes <- file.info(expected_panels)$size
  add_check(
    "All declared figure formats",
    all(file.exists(expected_panels)) &&
      all(is.finite(panel_sizes)) &&
      all(panel_sizes > 1000),
    paste(sum(file.exists(expected_panels)), "of", length(expected_panels))
  )
}

score_audit_path <- file.path(
  output_dir,
  "seven_donor_scores",
  "seven_donor_score_audit.csv"
)
if (file.exists(score_audit_path)) {
  score_audit <- read.csv(score_audit_path, check.names = FALSE)
  joined_cells <- score_audit$value[
    score_audit$metric == "Joined SFE cells"
  ]
  add_check(
    "Seven-donor score audit",
    length(joined_cells) == 1L &&
      joined_cells == 21271L &&
      all(score_audit$value > 0),
    paste(score_audit$metric, score_audit$value, collapse = "; ")
  )
  if (exists("sfe")) {
    expected_fib <- sum(
      as.character(sfe$lv1_anno_final) == "fibroblast",
      na.rm = TRUE
    )
    expected_ec <- sum(
      as.character(sfe$lv1_anno_final) == "endothelial",
      na.rm = TRUE
    )
    observed_fib <- score_audit$value[
      score_audit$metric == "Annotated fibroblast cells scored"
    ]
    observed_ec <- score_audit$value[
      score_audit$metric == "Annotated endothelial cells scored"
    ]
    add_check(
      "Score caches match finalized Fib/EC annotation",
      identical(as.numeric(observed_fib), as.numeric(expected_fib)) &&
        identical(as.numeric(observed_ec), as.numeric(expected_ec)),
      paste("fibroblast", observed_fib, "; endothelial", observed_ec)
    )
  }
}

validation_table <- bind_rows(checks)
validation_csv <- file.path(report_dir, "REPORT_VALIDATION.csv")
validation_txt <- file.path(report_dir, "REPORT_VALIDATION.txt")
write.csv(validation_table, validation_csv, row.names = FALSE)
writeLines(
  c(
    paste("Report directory:", report_dir),
    paste("Validation time:", format(Sys.time(), tz = "UTC")),
    "",
    paste0(
      ifelse(validation_table$passed, "PASS", "FAIL"),
      " | ",
      validation_table$check,
      " | ",
      validation_table$detail
    )
  ),
  validation_txt
)

print(validation_table, n = Inf)
if (any(!validation_table$passed)) {
  stop("One or more report validation checks failed.")
}
cat("\nAll report validation checks passed.\n")
