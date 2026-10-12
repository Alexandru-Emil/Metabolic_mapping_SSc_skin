# Yi-Nan Li: ResolviPrepFromSfe_20260224.Rmd; portable input adapter.
# Source: https://github.com/ynanli/ssc-metabolic/blob/f1bd384c1cb8a34b21ff2e238feec5e6a1d17cbc/scripts/ResolviPrepFromSfe_20260224.Rmd
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript prepare_resolvi_input.R input_sfe.rds output.h5ad")
suppressPackageStartupMessages({
  library(SpatialFeatureExperiment)
  library(SummarizedExperiment)
  library(Seurat)
  library(anndataR)
})
sfe <- readRDS(args[1])
if (!"counts" %in% assayNames(sfe)) stop("Input SFE must contain raw counts.")
if (!"SampleId" %in% colnames(colData(sfe))) stop("Input SFE must contain SampleId.")
samples <- as.character(sfe$SampleId)
if (anyNA(samples) || any(!nzchar(samples)) || length(unique(samples)) != 6L)
  stop("This workflow requires the original six-donor SFE after cell QC.")
counts <- assay(sfe, "counts")
cell_ids <- sub("^.*_", "", colnames(sfe))
if (anyDuplicated(cell_ids)) stop("Xenium cell IDs must be unique across the six-donor export.")
coords <- spatialCoords(sfe)
if (is.null(coords) || ncol(coords) != 2L || any(!is.finite(coords)))
  stop("Input SFE must contain finite two-dimensional spatial coordinates.")
# Resolve coordinate rows explicitly; do not silently reorder cells by membership.
if (is.null(rownames(coords))) {
  if (nrow(coords) != ncol(sfe)) stop("Spatial coordinate count does not match cells.")
} else {
  index <- match(colnames(sfe), rownames(coords))
  if (anyNA(index)) index <- match(cell_ids, rownames(coords))
  if (anyNA(index) || anyDuplicated(index)) stop("Spatial coordinates cannot be aligned uniquely to cells.")
  coords <- coords[index, , drop = FALSE]
}
meta <- data.frame(SampleId = samples, original_sfe_cell = colnames(sfe), row.names = colnames(sfe))
obj <- CreateSeuratObject(counts, meta.data = meta)
colnames(obj) <- cell_ids
obj$orig.ident <- samples
rownames(coords) <- cell_ids
colnames(coords) <- c("coord_1", "coord_2")
obj[["X_spatial"]] <- CreateDimReducObject(embeddings = as.matrix(coords), key = "coord_", assay = "RNA")
dir.create(dirname(args[2]), recursive = TRUE, showWarnings = FALSE)
# Explicit mappings guarantee the counts layer and the RESOLVI 1.4.1 spatial key.
write_h5ad(obj, args[2], x_mapping = "counts", layers_mapping = c(counts = "counts"),
           obs_mapping = TRUE, var_mapping = FALSE,
           obsm_mapping = c(X_spatial = "X_spatial"),
           varm_mapping = FALSE, obsp_mapping = FALSE, varp_mapping = FALSE, uns_mapping = FALSE)
