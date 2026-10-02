#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (requireNamespace("SpatialFeatureExperiment", quietly = TRUE)) {
    library(SpatialFeatureExperiment)
  }
  library(SpatialExperiment)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(zellkonverter)
})

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript export_sfe_for_mintflow.R --sfe-rds xenium_sfe.rds --cell-type-col cell_type [options]\n\n",
    "Options:\n",
    "  --out-dir DIR          Directory for MintFlow-ready h5ad files [mintflow_h5ad]\n",
    "  --slice-col COLUMN     colData column identifying tissue section/slice [sample_id]\n",
    "  --batch-col COLUMN     colData column identifying batch [defaults to --slice-col]\n",
    "  --assay-name NAME      Assay to write into AnnData .X; must be raw counts [counts]\n",
    sep = ""
  )
}

parse_cli <- function(args) {
  parsed <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("--help", "-h")) {
      usage()
      quit(save = "no", status = 0)
    }
    if (!startsWith(key, "--") || i == length(args) || startsWith(args[[i + 1L]], "--")) {
      stop("Expected --key value argument pair near: ", key, call. = FALSE)
    }
    parsed[[sub("^--", "", key)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  parsed
}

opt <- modifyList(
  list(
    "sfe-rds" = NULL,
    "out-dir" = "mintflow_h5ad",
    "cell-type-col" = NULL,
    "slice-col" = "sample_id",
    "batch-col" = NA_character_,
    "assay-name" = "counts"
  ),
  parse_cli(commandArgs(trailingOnly = TRUE))
)
if (is.null(opt$`sfe-rds`) || is.null(opt$`cell-type-col`)) {
  usage()
  stop("Required: --sfe-rds and --cell-type-col", call. = FALSE)
}

batch_col <- if (is.na(opt$`batch-col`)) opt$`slice-col` else opt$`batch-col`
sfe <- readRDS(opt$`sfe-rds`)

if (!opt$`assay-name` %in% assayNames(sfe)) {
  stop("Assay not found: ", opt$`assay-name`, call. = FALSE)
}

cd <- as.data.frame(colData(sfe))
if (!opt$`cell-type-col` %in% names(cd)) {
  stop("Cell type column not found in colData: ", opt$`cell-type-col`, call. = FALSE)
}
if (!opt$`slice-col` %in% names(cd)) {
  cd[[opt$`slice-col`]] <- "section1"
}
if (!batch_col %in% names(cd)) {
  cd[[batch_col]] <- cd[[opt$`slice-col`]]
}

coords <- spatialCoords(sfe)
if (ncol(coords) < 2) {
  stop("spatialCoords(sfe) must contain at least two coordinate columns.", call. = FALSE)
}

cd$x_centroid <- as.numeric(coords[, 1])
cd$y_centroid <- as.numeric(coords[, 2])
cd$celltype_for_MintFlow <- as.factor(cd[[opt$`cell-type-col`]])
cd$TissueSectionID_for_MintFlow <- as.factor(cd[[opt$`slice-col`]])
cd$batchID_for_MintFlow <- as.factor(cd[[batch_col]])
colData(sfe) <- DataFrame(cd, row.names = colnames(sfe))

dir.create(opt$`out-dir`, recursive = TRUE, showWarnings = FALSE)
indices_by_section <- split(seq_len(ncol(sfe)), cd$TissueSectionID_for_MintFlow)
out_files <- character(length(indices_by_section))

for (i in seq_along(indices_by_section)) {
  sce <- as(sfe[, indices_by_section[[i]]], "SingleCellExperiment")
  out_files[i] <- file.path(opt$`out-dir`, sprintf("tissue_section_%03d.h5ad", i))
  writeH5AD(sce, out_files[i], X_name = opt$`assay-name`)
}

out_files <- normalizePath(out_files, winslash = "/", mustWork = FALSE)
writeLines(out_files, file.path(opt$`out-dir`, "h5ad_files.txt"))
message("Wrote ", length(out_files), " h5ad file(s) and h5ad_files.txt to ", opt$`out-dir`)
