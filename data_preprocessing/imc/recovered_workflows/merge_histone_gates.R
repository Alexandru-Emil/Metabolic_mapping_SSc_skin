# Recovered from histone merged.R (second, post-Xenium IMC cohort).
merge_histone_gates <- function(spe, gate_objects, label_column = "gates") {
  ids <- colnames(spe)
  if (is.null(ids) || anyDuplicated(ids)) stop("Full-object cell IDs must be unique.")
  cd <- SummarizedExperiment::colData(spe)
  if (!label_column %in% colnames(cd)) stop("Missing label column.")
  positive_ids <- unique(unlist(lapply(gate_objects, colnames), use.names = FALSE))
  if (!length(positive_ids)) stop("No Histone-positive cells supplied.")
  if (any(!positive_ids %in% ids)) stop("Histone gates contain unknown cell IDs.")
  labels <- as.character(cd[[label_column]])
  if (anyNA(labels)) stop("Missing input labels.")
  labels[labels == "unlabeled" & ids %in% positive_ids] <- "Stromal"
  cd[[label_column]] <- labels
  SummarizedExperiment::colData(spe) <- cd
  spe
}
if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  if (length(args) != 3L) stop("Usage: Rscript merge_histone_gates.R INPUT.rds GATE_DIR OUTPUT.rds")
  suppressPackageStartupMessages(library(SpatialExperiment))
  files <- list.files(args[2L], pattern = "\\.rds$", full.names = TRUE,
                      recursive = TRUE, ignore.case = TRUE)
  if (!length(files)) stop("No RDS gate exports found.")
  saveRDS(merge_histone_gates(readRDS(args[1L]), lapply(files, readRDS)), args[3L])
}
