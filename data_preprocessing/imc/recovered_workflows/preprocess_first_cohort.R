# Recovered from the 2023 first-cohort Steinbock import workflow.
# Input: an unfiltered SpatialExperiment containing the counts assay and area.
preprocess_first_cohort <- function(spe, lower = 20, upper = 250, cofactor = 1) {
  stopifnot(inherits(spe, "SpatialExperiment"), cofactor > 0, lower <= upper)
  area <- SummarizedExperiment::colData(spe)$area
  if (is.null(area) || anyNA(area)) stop("A complete area column is required.")
  keep <- area >= lower & area <= upper
  SummarizedExperiment::colData(spe)$size.filter <- ifelse(keep, "include", "exclude")
  spe <- spe[, keep]
  x <- asinh(SummarizedExperiment::assay(spe, "counts") / cofactor)
  SummarizedExperiment::assay(spe, "asinh") <- x
  SummarizedExperiment::assay(spe, "rescaled") <- t(apply(x, 1L, function(v) {
    (v - min(v)) / (max(v) - min(v))
  }))
  SummarizedExperiment::assay(spe, "z-rescaled") <- t(apply(x, 1L, function(v) {
    (v - mean(v)) / sd(v)
  }))
  spe
}
if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  if (length(args) != 2L) stop("Usage: Rscript preprocess_first_cohort.R INPUT.rds OUTPUT.rds")
  suppressPackageStartupMessages(library(SpatialExperiment))
  saveRDS(preprocess_first_cohort(readRDS(args[1L])), args[2L])
}
