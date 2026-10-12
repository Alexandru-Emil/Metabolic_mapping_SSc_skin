# Recovered from Script to merge and subset CD45 gate.txt.
# Only CD45-positive cells receive the four immune sublabels.
merge_cd68_cd3_gates <- function(spe, gate_objects) {
  ids <- colnames(spe)
  if (is.null(ids) || anyDuplicated(ids)) stop("Full-object cell IDs must be unique.")
  cd <- SummarizedExperiment::colData(spe)
  labels <- as.character(cd$cell_labels)
  if (length(labels) != ncol(spe) || anyNA(labels)) stop("Complete cell_labels required.")
  positive <- list(CD68 = character(), CD3 = character())
  for (g in gate_objects) {
    lab <- as.character(SummarizedExperiment::colData(g)$cytomapper_CellLabel)
    if (length(lab) != ncol(g) || anyNA(lab)) stop("Invalid gate labels.")
    if (any(!lab %in% names(positive))) stop("Expected CD68 or CD3 gate labels.")
    if (any(!colnames(g) %in% ids)) stop("Gate contains unknown cell IDs.")
    for (name in names(positive)) positive[[name]] <- union(positive[[name]], colnames(g)[lab == name])
  }
  if (!length(gate_objects)) stop("No gate objects supplied.")
  p68 <- ids %in% positive$CD68
  p3 <- ids %in% positive$CD3
  immune <- labels == "CD45+"
  labels[immune & !p68 & !p3] <- "CD45+CD68-CD3-"
  labels[immune & p68 & !p3] <- "CD45+CD68+"
  labels[immune & !p68 & p3] <- "CD45+CD3+"
  labels[immune & p68 & p3] <- "CD45+CD68+CD3+"
  cd$new_celllabels <- labels
  SummarizedExperiment::colData(spe) <- cd
  spe
}
if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  if (length(args) != 3L) stop("Usage: Rscript merge_cd68_cd3_gates.R INPUT.rds GATE_DIR OUTPUT.rds")
  suppressPackageStartupMessages(library(SpatialExperiment))
  files <- list.files(args[2L], pattern = "\\.rds$", full.names = TRUE, ignore.case = TRUE)
  if (!length(files)) stop("No RDS gate exports found.")
  saveRDS(merge_cd68_cd3_gates(readRDS(args[1L]), lapply(files, readRDS)), args[3L])
}
