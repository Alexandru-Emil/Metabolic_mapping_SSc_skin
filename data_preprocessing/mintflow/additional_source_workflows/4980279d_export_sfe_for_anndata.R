#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SpatialFeatureExperiment)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: export_sfe_for_anndata.R <sfe.rds> <spe.RData> <out_dir>")
}

path_sfe <- args[[1]]
path_spe <- args[[2]]
out_dir <- args[[3]]

to_dgC <- function(x) {
  if (inherits(x, "dgCMatrix")) {
    return(x)
  }
  as(as(x, "CsparseMatrix"), "dgCMatrix")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sfe <- readRDS(path_sfe)
env <- new.env(parent = emptyenv())
load(path_spe, envir = env)
if (!exists("spe", envir = env)) {
  stop("Expected object named 'spe' in ", path_spe)
}
spe <- env$spe

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
sfe$x_centroid <- as.numeric(coords[, "x_centroid"])
sfe$y_centroid <- as.numeric(coords[, "y_centroid"])

counts <- to_dgC(assay(sfe, "counts"))
storage.mode(counts@x) <- "double"

metadata <- as.data.frame(colData(sfe))
metadata$cell_id <- colnames(sfe)
metadata <- metadata[, c("cell_id", setdiff(colnames(metadata), "cell_id")), drop = FALSE]

genes <- data.frame(gene_id = rownames(sfe), gene_symbol = rownames(sfe), stringsAsFactors = FALSE)

Matrix::writeMM(counts, file.path(out_dir, "counts_genes_by_cells.mtx"))
write.csv(metadata, file.path(out_dir, "obs.csv"), row.names = FALSE, quote = TRUE)
write.csv(genes, file.path(out_dir, "var.csv"), row.names = FALSE, quote = TRUE)

message("Exported Matrix Market bundle to ", out_dir)
message("Cells: ", ncol(sfe), " Genes: ", nrow(sfe))
print(table(sfe$metfiblabel, useNA = "ifany"))
