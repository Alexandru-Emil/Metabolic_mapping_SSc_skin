# Original matching convention: Yi-Nan Li, 20260318_ImcMethiFib.Rmd.
# Attach the six-donor Seurat annotations to their input SFE without positional joins.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript attach_annotations.R input_sfe.rds annotated_seurat.rds output_sfe.rds")
suppressPackageStartupMessages({library(SpatialFeatureExperiment); library(Seurat)})
sfe <- readRDS(args[1])
obj <- readRDS(args[2])
meta <- obj[[]]
if (!"SampleId" %in% colnames(meta) || anyNA(meta$SampleId) || length(unique(meta$SampleId)) != 6L)
  stop("The annotated object must contain the original six donors.")
if ("original_sfe_cell" %in% colnames(meta)) {
  keys <- as.character(meta$original_sfe_cell)
} else {
  # Historical model outputs use bare Xenium barcodes; orig.ident is the sample.
  if (!"orig.ident" %in% colnames(meta)) stop("Missing original_sfe_cell or orig.ident mapping.")
  keys <- paste0(meta$orig.ident, "_", colnames(obj))
}
if (anyNA(keys) || anyDuplicated(keys)) stop("Annotation keys must be complete and unique.")
index <- match(colnames(sfe), keys)
if (anyNA(index)) stop("Use the same QC-filtered SFE used to prepare the model input; some cells do not match.")
if (!all(as.character(sfe$SampleId) == as.character(meta$SampleId[index])))
  stop("Sample labels disagree between the SFE and annotation object.")
for (field in c("lv1_anno", "lv2_anno")) {
  if (!field %in% colnames(meta)) stop(paste("Missing annotation field", field))
  SummarizedExperiment::colData(sfe)[[field]] <- as.character(meta[[field]][index])
}
for (reduction in intersect(c("X_resolvi", "umap"), Reductions(obj))) {
  embedding <- Embeddings(obj, reduction = reduction)
  SingleCellExperiment::reducedDim(sfe, reduction) <- embedding[index, , drop = FALSE]
}
dir.create(dirname(args[3]), recursive = TRUE, showWarnings = FALSE)
saveRDS(sfe, args[3])
