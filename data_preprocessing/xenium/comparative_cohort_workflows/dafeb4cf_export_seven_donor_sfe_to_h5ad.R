#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SpatialFeatureExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(Seurat)
  library(anndataR)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    "Usage: export_seven_donor_sfe_to_h5ad.R ",
    "<input_sfe.rds> <output.h5ad>"
  )
}

input_path <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_path <- normalizePath(
  args[[2]],
  winslash = "/",
  mustWork = FALSE
)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

donor_levels <- c(
  "Validation1",
  "Validation2",
  "ExcludedValidation",
  "Validation3",
  "Validation4",
  "Validation5",
  "Validation6"
)

# These acquisition-slide assignments come from the original ResolVI input
# H5AD created by 20251118_ResolviPrep.Rmd.
slide_by_donor <- c(
  "Validation1" = "0037406",
  "Validation2" = "0037406",
  "ExcludedValidation" = "0037085",
  "Validation3" = "0037406",
  "Validation4" = "0037406",
  "Validation5" = "0037085",
  "Validation6" = "0037085"
)

message("Reading seven-donor SFE: ", input_path)
sfe <- readRDS(input_path)
stopifnot(
  inherits(sfe, "SpatialFeatureExperiment"),
  identical(sort(unique(as.character(sfe$SampleId))), sort(donor_levels)),
  nrow(sfe) == 5099L,
  ncol(sfe) == 21271L,
  !anyDuplicated(rownames(sfe)),
  !anyDuplicated(colnames(sfe)),
  "counts" %in% assayNames(sfe)
)

counts <- assay(sfe, "counts")
if (!inherits(counts, "dgCMatrix")) {
  counts <- as(as(counts, "CsparseMatrix"), "dgCMatrix")
}
stopifnot(
  Matrix::nnzero(counts) > 0,
  all(counts@x >= 0),
  all(counts@x == round(counts@x)),
  all(Matrix::colSums(counts) > 0)
)

coords <- as.matrix(spatialCoords(sfe))
coords <- coords[, seq_len(2L), drop = FALSE]
storage.mode(coords) <- "double"
rownames(coords) <- colnames(sfe)
colnames(coords) <- c("coord_1", "coord_2")
stopifnot(
  nrow(coords) == ncol(sfe),
  all(is.finite(coords)),
  identical(rownames(coords), colnames(sfe))
)

source_metadata <- as.data.frame(colData(sfe), check.names = FALSE)
metadata_columns <- intersect(
  c(
    "SampleId",
    "sample_id",
    "ROI",
    "Patient",
    "Group",
    "KEY",
    "imc_matched",
    "metfiblabel",
    "metEClabel",
    "lv1_anno_new",
    "lv2_anno_new",
    "lv1_anno_joint",
    "lv2_anno_joint",
    "lv1_anno_mapped",
    "lv2_anno_mapped",
    "lv1_mapping_score",
    "lv2_mapping_score",
    "lv1_anno_final",
    "lv2_anno_final",
    "annotation_source",
    "annotation_strategy_final"
  ),
  colnames(source_metadata)
)
obs <- source_metadata[, metadata_columns, drop = FALSE]
for (column in colnames(obs)) {
  if (is.factor(obs[[column]])) {
    obs[[column]] <- as.character(obs[[column]])
  }
}
obs$SampleId <- as.character(obs$SampleId)
obs$SlideId <- unname(slide_by_donor[obs$SampleId])
obs$cell_id <- colnames(sfe)
rownames(obs) <- colnames(sfe)

stopifnot(
  !anyNA(obs$SlideId),
  identical(names(table(obs$SampleId)), donor_levels),
  identical(rownames(obs), colnames(sfe))
)

message("Building the Seurat export object")
seurat_object <- CreateSeuratObject(
  counts = counts,
  meta.data = obs,
  project = "seven_donor_joint_resolvi"
)
seurat_object[["X_spatial"]] <- CreateDimReducObject(
  embeddings = coords,
  key = "coord_",
  assay = "RNA"
)

message("Writing H5AD: ", output_path)
if (file.exists(output_path)) {
  file.remove(output_path)
}
write_h5ad(
  seurat_object,
  output_path,
  compression = "gzip"
)

cat("\nExport complete\n")
cat("Input:", input_path, "\n")
cat("Output:", normalizePath(output_path, winslash = "/"), "\n")
cat("Dimensions:", ncol(sfe), "cells x", nrow(sfe), "genes\n")
cat("Nonzero counts:", Matrix::nnzero(counts), "\n")
print(table(obs$SampleId, obs$SlideId))
