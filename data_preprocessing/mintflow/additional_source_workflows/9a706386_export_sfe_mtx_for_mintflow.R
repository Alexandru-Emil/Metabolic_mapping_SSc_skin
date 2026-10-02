#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (requireNamespace("SpatialFeatureExperiment", quietly = TRUE)) {
    library(SpatialFeatureExperiment)
  }
  library(SpatialExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(Matrix)
})

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript export_sfe_mtx_for_mintflow.R --sfe-rds xenium_sfe.rds --cell-type-col cell_type [options]\n\n",
    "Options:\n",
    "  --out-dir DIR          Directory for Matrix Market exports [mintflow_mtx]\n",
    "  --slice-col COLUMN     colData column identifying tissue section/slice [sample_id]\n",
    "  --batch-col COLUMN     colData column identifying batch [defaults to --slice-col]\n",
    "  --assay-name NAME      Assay to export; must be raw counts [counts]\n",
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
    "out-dir" = "mintflow_mtx",
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
cd$celltype_for_MintFlow <- as.character(cd[[opt$`cell-type-col`]])
cd$TissueSectionID_for_MintFlow <- as.character(cd[[opt$`slice-col`]])
cd$batchID_for_MintFlow <- as.character(cd[[batch_col]])

if (anyNA(cd$celltype_for_MintFlow)) {
  stop("Cell type column contains NA values; choose or fill a complete annotation column.", call. = FALSE)
}

counts <- assay(sfe, opt$`assay-name`)
if (!inherits(counts, "sparseMatrix")) {
  counts <- as(counts, "dgCMatrix")
}

feature_ids <- rownames(sfe)
if (is.null(feature_ids)) {
  feature_ids <- sprintf("feature_%06d", seq_len(nrow(sfe)))
}
cell_ids <- colnames(sfe)
if (is.null(cell_ids)) {
  cell_ids <- sprintf("cell_%06d", seq_len(ncol(sfe)))
}
if (anyDuplicated(cell_ids)) {
  stop("Cell IDs/colnames are duplicated; AnnData obs names must be unique.", call. = FALSE)
}

dir.create(opt$`out-dir`, recursive = TRUE, showWarnings = FALSE)
indices_by_section <- split(seq_len(ncol(sfe)), cd$TissueSectionID_for_MintFlow)
manifest <- data.frame(section = character(), directory = character(), stringsAsFactors = FALSE)

for (i in seq_along(indices_by_section)) {
  section_dir <- file.path(opt$`out-dir`, sprintf("tissue_section_%03d", i))
  dir.create(section_dir, recursive = TRUE, showWarnings = FALSE)

  idx <- indices_by_section[[i]]
  Matrix::writeMM(t(counts[, idx, drop = FALSE]), file.path(section_dir, "counts_cells_by_genes.mtx"))

  obs <- cd[idx, , drop = FALSE]
  rownames(obs) <- cell_ids[idx]
  write.csv(obs, file.path(section_dir, "obs.csv"), quote = TRUE)

  var <- as.data.frame(rowData(sfe))
  if (ncol(var) == 0L) {
    var <- data.frame(gene_id = feature_ids, row.names = feature_ids)
  } else {
    rownames(var) <- feature_ids
    var$gene_id <- feature_ids
  }
  write.csv(var, file.path(section_dir, "var.csv"), quote = TRUE)

  manifest <- rbind(
    manifest,
    data.frame(section = names(indices_by_section)[[i]], directory = normalizePath(section_dir, winslash = "/", mustWork = FALSE))
  )
}

write.csv(manifest, file.path(opt$`out-dir`, "manifest.csv"), row.names = FALSE, quote = TRUE)
message("Wrote ", nrow(manifest), " Matrix Market section export(s) to ", opt$`out-dir`)
