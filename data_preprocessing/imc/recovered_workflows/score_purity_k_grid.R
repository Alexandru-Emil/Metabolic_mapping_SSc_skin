# Recover purity/silhouette scoring and selection from the 2023-12-05 history.
# Clustering itself is deliberately an input: saved Rphenograph_cluster_<k> columns.
score_purity_k_grid <- function(spe, markers, k_grid = c(5,20,50,100,150,200,250,500),
                                assay_name = "z_c_rescaled", label_column = "cell_labels") {
  x <- SummarizedExperiment::assay(spe, assay_name)
  cd <- SummarizedExperiment::colData(spe)
  if (is.character(markers)) markers <- match(markers, rownames(spe))
  if (!length(markers) || anyNA(markers) || anyDuplicated(markers) ||
      any(markers < 1 | markers > nrow(spe))) stop("Invalid marker selection.")
  types <- as.character(cd[[label_column]])
  if (length(types) != ncol(spe) || anyNA(types)) stop("Complete compartment labels required.")
  out <- list()
  for (type in unique(types)) {
    if (type == "unlabeled") next
    selected <- which(types == type)
    mat <- t(x[markers, selected, drop = FALSE])
    if (any(!is.finite(mat))) stop("Nonfinite scoring values.")
    for (k in k_grid) {
      method <- paste0("Rphenograph_cluster_", k)
      cl <- cd[[method]][selected]
      if (length(cl) != length(selected) || anyNA(cl)) stop("Missing clustering: ", method)
      # Historical defaults: neighborhood k=50; centroid-based approximate silhouette.
      p <- bluster::neighborPurity(mat, cl, k = 50)$purity
      s <- bluster::approxSilhouette(mat, cl)$width
      pc <- vapply(split(p, cl), mean, numeric(1))
      sc <- vapply(split(s, cl), mean, numeric(1))
      out[[length(out) + 1L]] <- data.frame(
        cell_type = type, cluster_method = method, n_clusters = length(pc),
        mean_purity_cells = mean(p), mean_purity_cluster_means = mean(pc),
        median_purity_cluster_means = median(pc), min_purity_cluster_means = min(pc),
        max_purity_cluster_means = max(pc), mean_silhouette_cells = mean(s),
        mean_silhouette_cluster_means = mean(sc), median_silhouette_cluster_means = median(sc),
        min_silhouette_cluster_means = min(sc), max_silhouette_cluster_means = max(sc))
    }
  }
  do.call(rbind, out)
}
legacy_purity_table <- function(scores) {
  # Preserve the archive's column order, mislabeled headings, and 3-decimal rounding.
  names(scores) <- c("Cell type", "Cluster method", "number of cluster",
    "mean purity across cells", "median purity across cells", "mean purity across clusters",
    "minimum purity across clusters", "maximum purity across clusters",
    "mean silhouette across cells", "median silhouette across cells", "mean silhouette across clusters",
    "minimum silhouette across clusters", "maximum silhouette across clusters")
  scores[4:13] <- lapply(scores[4:13], round, digits = 3)
  scores
}
select_purity_k <- function(scores, historical_rounding = TRUE) {
  if (historical_rounding) scores[4:13] <- lapply(scores[4:13], round, digits = 3)
  rows <- lapply(unique(scores$cell_type), function(type) {
    s <- scores[scores$cell_type == type, , drop = FALSE]
    total <- s$mean_purity_cells + s$mean_silhouette_cells
    best <- which(total == max(total))
    if (length(best) != 1L) stop("Tied optimum for ", type, "; inspect all candidates.")
    data.frame(cell_type = type, cluster_method = s$cluster_method[best])
  })
  do.call(rbind, rows)
}
if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  if (length(args) != 3L) stop("Usage: Rscript score_purity_k_grid.R INPUT.rds MARKERS.txt OUTPUT_PREFIX")
  suppressPackageStartupMessages(library(SpatialExperiment))
  markers <- readLines(args[2L], warn = FALSE)
  scores <- score_purity_k_grid(readRDS(args[1L]), markers)
  write.csv(scores, paste0(args[3L], "_scores.csv"), row.names = FALSE)
  write.csv(legacy_purity_table(scores), paste0(args[3L], "_legacy.csv"), row.names = FALSE)
  write.csv(select_purity_k(scores), paste0(args[3L], "_selected.csv"), row.names = FALSE)
}
