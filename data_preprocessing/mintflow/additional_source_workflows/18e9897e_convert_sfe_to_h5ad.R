#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SpatialFeatureExperiment)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(zellkonverter)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: convert_sfe_to_h5ad.R <sfe.rds> <spe.RData> <out.h5ad>")
}

path_sfe <- args[[1]]
path_spe <- args[[2]]
path_out <- args[[3]]

to_dgC <- function(x) {
  if (inherits(x, "dgCMatrix")) {
    return(x)
  }
  as(as(x, "CsparseMatrix"), "dgCMatrix")
}

sfe <- readRDS(path_sfe)
env <- new.env(parent = emptyenv())
load(path_spe, envir = env)
if (!exists("spe", envir = env)) {
  stop("Expected object named 'spe' in ", path_spe)
}
spe <- env$spe

if (!"KEY" %in% colnames(colData(sfe))) {
  stop("The SFE object does not contain colData(sfe)$KEY for IMC label transfer.")
}
if (!"Rphenograph_500" %in% colnames(colData(spe))) {
  stop("The IMC SPE object does not contain colData(spe)$Rphenograph_500.")
}

spe_rphenograph <- as.character(colData(spe)[["Rphenograph_500"]])
metfiblabel <- ifelse(
  spe_rphenograph %in% c("1", "2"),
  "Met_hi_Fib",
  ifelse(spe_rphenograph == "3", "Other_Fib", NA_character_)
)

idx <- match(colData(sfe)[["KEY"]], colnames(spe))
sfe$metfiblabel <- metfiblabel[idx]

keep <- !is.na(sfe$metfiblabel) & sfe$metfiblabel != ""
sfe <- sfe[, keep]
sfe$metfiblabel <- factor(as.character(sfe$metfiblabel))
sfe$mintflow_slice_id <- "xenium_sfe_09032026"
sfe$mintflow_batch <- "xenium_sfe_09032026"

coords <- spatialCoords(sfe)
if (!all(c("x_centroid", "y_centroid") %in% colnames(coords))) {
  stop("Expected spatial coordinates named x_centroid and y_centroid.")
}
sfe$x_centroid <- as.numeric(coords[, "x_centroid"])
sfe$y_centroid <- as.numeric(coords[, "y_centroid"])

counts <- to_dgC(assay(sfe, "counts"))
storage.mode(counts@x) <- "integer"

sce <- SingleCellExperiment(
  assays = list(counts = counts),
  rowData = rowData(sfe),
  colData = colData(sfe)
)
reducedDim(sce, "spatial") <- as.matrix(coords[, c("x_centroid", "y_centroid")])

dir.create(dirname(path_out), recursive = TRUE, showWarnings = FALSE)
writeH5AD(sce, path_out, X_name = "counts")

message("Wrote ", path_out)
message("Cells: ", ncol(sce), " Genes: ", nrow(sce))
message("metfiblabel:")
print(table(sce$metfiblabel, useNA = "ifany"))
