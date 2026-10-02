source("data_preprocessing/common/config.R")
suppressPackageStartupMessages({
  library(SpatialFeatureExperiment)
  library(SpatialExperiment)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(S4Vectors)
  library(Matrix)
  library(scuttle)
  library(dplyr)
})

set.seed(20260729)

project_root <- normalizePath(
  file.path(getwd(), "seven_donor_reanalysis"),
  winslash = "/",
  mustWork = FALSE
)
output_dir <- file.path(project_root, "outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sfe6_path <- paste0(
  paste0(project_path("Revision/"), "/"),
  "Nan annotation/Data/Xenium sfe objects/Xenium sfe object/sfe_09032026.rds"
)
sfe230_path <- paste0(
  paste0(project_path("Revision/"), "/"),
  "Nan annotation/Data/Xenium sfe objects/XeniumSfe_ExcludedValidation_260729115149.rds"
)
spe_path <- paste0(
  paste0(project_path("Revision/"), "/"),
  "Nan annotation/Data/IMC spe objects/New IMC/spe2025-12-17 umap.RData"
)

joined_path <- file.path(output_dir, "sfe_seven_donors_joined_unannotated.rds")

to_dgC <- function(x) {
  if (inherits(x, "dgCMatrix")) {
    return(x)
  }
  as(as(x, "CsparseMatrix"), "dgCMatrix")
}

common_type <- function(x, y) {
  classes <- c(class(x), class(y))
  if (any(classes %in% c("character", "factor"))) return("character")
  if (any(classes %in% c("numeric", "integer"))) return("numeric")
  if (all(classes %in% "logical")) return("logical")
  "character"
}

cast_column <- function(x, type) {
  switch(
    type,
    character = as.character(x),
    numeric = as.numeric(x),
    logical = as.logical(x),
    as.character(x)
  )
}

missing_column <- function(n, type) {
  switch(
    type,
    character = rep(NA_character_, n),
    numeric = rep(NA_real_, n),
    logical = rep(NA, n),
    rep(NA_character_, n)
  )
}

harmonize_coldata <- function(lhs, rhs) {
  lhs_df <- as.data.frame(colData(lhs), check.names = FALSE)
  rhs_df <- as.data.frame(colData(rhs), check.names = FALSE)
  all_names <- union(colnames(lhs_df), colnames(rhs_df))

  lhs_out <- vector("list", length(all_names))
  rhs_out <- vector("list", length(all_names))
  names(lhs_out) <- all_names
  names(rhs_out) <- all_names

  for (nm in all_names) {
    lhs_has <- nm %in% colnames(lhs_df)
    rhs_has <- nm %in% colnames(rhs_df)

    if (lhs_has && rhs_has) {
      type <- common_type(lhs_df[[nm]], rhs_df[[nm]])
    } else if (lhs_has) {
      type <- if (is.factor(lhs_df[[nm]])) "character" else
        if (is.integer(lhs_df[[nm]]) || is.numeric(lhs_df[[nm]])) "numeric" else
          if (is.logical(lhs_df[[nm]])) "logical" else "character"
    } else {
      type <- if (is.factor(rhs_df[[nm]])) "character" else
        if (is.integer(rhs_df[[nm]]) || is.numeric(rhs_df[[nm]])) "numeric" else
          if (is.logical(rhs_df[[nm]])) "logical" else "character"
    }

    lhs_out[[nm]] <- if (lhs_has) {
      cast_column(lhs_df[[nm]], type)
    } else {
      missing_column(nrow(lhs_df), type)
    }
    rhs_out[[nm]] <- if (rhs_has) {
      cast_column(rhs_df[[nm]], type)
    } else {
      missing_column(nrow(rhs_df), type)
    }
  }

  lhs_out <- DataFrame(lhs_out, row.names = rownames(lhs_df), check.names = FALSE)
  rhs_out <- DataFrame(rhs_out, row.names = rownames(rhs_df), check.names = FALSE)
  list(lhs = lhs_out, rhs = rhs_out)
}

combine_sf_geometry <- function(lhs, rhs) {
  lhs_df <- lhs
  rhs_df <- rhs
  geom_lhs <- attr(lhs_df, "sf_column")
  geom_rhs <- attr(rhs_df, "sf_column")

  if (!identical(geom_lhs, geom_rhs)) {
    names(rhs_df)[names(rhs_df) == geom_rhs] <- geom_lhs
    sf::st_geometry(rhs_df) <- geom_lhs
  }

  all_names <- union(colnames(lhs_df), colnames(rhs_df))
  for (nm in setdiff(all_names, colnames(lhs_df))) {
    lhs_df[[nm]] <- NA
  }
  for (nm in setdiff(all_names, colnames(rhs_df))) {
    rhs_df[[nm]] <- NA
  }

  lhs_df <- lhs_df[, all_names, drop = FALSE]
  rhs_df <- rhs_df[, all_names, drop = FALSE]
  out <- rbind(lhs_df, rhs_df)
  rownames(out) <- c(rownames(lhs), rownames(rhs))
  out
}

message("Reading six-donor SFE")
sfe6 <- readRDS(sfe6_path)
message("Reading portable ExcludedValidation SFE")
sfe230_all <- readRDS(sfe230_path)

stopifnot(
  inherits(sfe6, "SpatialFeatureExperiment"),
  inherits(sfe230_all, "SpatialFeatureExperiment"),
  identical(rownames(sfe6), rownames(sfe230_all)),
  identical(as.data.frame(rowData(sfe6)), as.data.frame(rowData(sfe230_all))),
  !anyDuplicated(colnames(sfe6)),
  !anyDuplicated(colnames(sfe230_all)),
  length(intersect(colnames(sfe6), colnames(sfe230_all))) == 0L
)

counts230 <- assay(sfe230_all, "counts")
stopifnot(
  inherits(counts230, "dgCMatrix"),
  Matrix::nnzero(counts230) > 0,
  sum(counts230) > 0
)

qc230 <- as.logical(colData(sfe230_all)[["qc_pass"]])
qc230[is.na(qc230)] <- FALSE
sfe230 <- sfe230_all[, qc230]

load(spe_path)
stopifnot(exists("spe"), "KEY" %in% colnames(colData(sfe230_all)))

key230_all <- as.character(colData(sfe230_all)[["KEY"]])
key230_keep <- as.character(colData(sfe230)[["KEY"]])

qc_audit <- data.frame(
  metric = c(
    "ExcludedValidation cells supplied",
    "ExcludedValidation qc_pass cells retained",
    "ExcludedValidation qc_fail cells excluded",
    "ExcludedValidation cells matched to IMC before QC",
    "ExcludedValidation cells matched to IMC after QC"
  ),
  n_cells = c(
    ncol(sfe230_all),
    ncol(sfe230),
    sum(!qc230),
    sum(key230_all %in% colnames(spe)),
    sum(key230_keep %in% colnames(spe))
  ),
  stringsAsFactors = FALSE
)
write.csv(qc_audit, file.path(output_dir, "ssc230_qc_join_audit.csv"), row.names = FALSE)
print(qc_audit)

message("Reducing both inputs to the portable assays needed downstream")
assays(sfe6) <- SimpleList(counts = to_dgC(assay(sfe6, "counts")))
assays(sfe230) <- SimpleList(counts = to_dgC(assay(sfe230, "counts")))

reducedDims(sfe6) <- SimpleList()
reducedDims(sfe230) <- SimpleList()

cd_pair <- harmonize_coldata(sfe6, sfe230)
colData(sfe6) <- cd_pair$lhs
colData(sfe230) <- cd_pair$rhs

message("Combining assays, metadata, and sf geometries through the public constructor")
counts7 <- cbind(
  assay(sfe6, "counts"),
  assay(sfe230, "counts")
)
coldata7 <- rbind(colData(sfe6), colData(sfe230))

centroids6 <- colGeometries(sfe6)[["centroids"]]
centroids230 <- colGeometries(sfe230)[["centroids"]]
cellseg6 <- colGeometries(sfe6)[["cellSeg"]]
cellseg230 <- colGeometries(sfe230)[["cellSeg"]]

rownames(centroids6) <- colnames(sfe6)
rownames(centroids230) <- colnames(sfe230)
rownames(cellseg6) <- colnames(sfe6)
rownames(cellseg230) <- colnames(sfe230)

centroids7 <- combine_sf_geometry(centroids6, centroids230)
cellseg7 <- combine_sf_geometry(cellseg6, cellseg230)

message(
  "Combined dimensions: counts=", ncol(counts7),
  ", colData=", nrow(coldata7),
  ", centroids=", nrow(centroids7),
  ", cellSeg=", nrow(cellseg7)
)
stopifnot(
  ncol(counts7) == nrow(coldata7),
  ncol(counts7) == nrow(centroids7),
  ncol(counts7) == nrow(cellseg7)
)

sfe7 <- SpatialFeatureExperiment(
  assays = list(counts = counts7),
  rowData = rowData(sfe6),
  colData = coldata7,
  sample_id = as.character(coldata7[["SampleId"]]),
  spatialCoordsNames = c("x_centroid", "y_centroid"),
  colGeometries = list(
    centroids = centroids7,
    cellSeg = cellseg7
  ),
  unit = "micron"
)

stopifnot(
  nrow(sfe7) == 5099L,
  ncol(sfe7) == ncol(sfe6) + ncol(sfe230),
  length(unique(as.character(sfe7$SampleId))) == 7L,
  !anyDuplicated(colnames(sfe7)),
  all(c("centroids", "cellSeg") %in% colGeometryNames(sfe7))
)

message("Recalculating size factors and logcounts jointly across seven donors")
SingleCellExperiment::sizeFactors(sfe7) <- NULL
sfe7 <- scuttle::logNormCounts(sfe7)
assay(sfe7, "counts") <- to_dgC(assay(sfe7, "counts"))
assay(sfe7, "logcounts") <- to_dgC(assay(sfe7, "logcounts"))

sfe7$cohort_donor <- as.character(sfe7$SampleId)
sfe7$cohort_donor <- factor(
  sfe7$cohort_donor,
  levels = c("Validation1", "Validation2", "ExcludedValidation", "Validation3", "Validation4", "Validation5", "Validation6")
)
sfe7$annotation_source <- ifelse(
  as.character(sfe7$SampleId) == "ExcludedValidation",
  "ExcludedValidation pending seven-donor annotation",
  "Finalized six-donor ResolVI annotation"
)

metadata(sfe7)$seven_donor_join <- list(
  created = as.character(Sys.time()),
  six_donor_source = sfe6_path,
  ssc230_source = sfe230_path,
  imc_source = spe_path,
  ssc230_qc_rule = "Retained qc_pass == TRUE to match the existing six-donor SFE",
  assays_retained = c("counts", "logcounts"),
  assays_omitted = c(
    "H0",
    "H1"
  ),
  omission_reason = paste(
    "H0/H1 were legacy spatial helper assays absent from ExcludedValidation and are not",
    "used by the clustering, annotation, scMetabolism, ECM, or figure workflows."
  )
)

donor_audit <- as.data.frame(colData(sfe7)) %>%
  transmute(
    donor = as.character(SampleId),
    KEY = as.character(KEY),
    qc_pass = as.logical(qc_pass)
  ) %>%
  mutate(
    valid_key = !is.na(KEY) & nzchar(KEY) & KEY != "not found_not found",
    matched_imc = !is.na(KEY) & KEY %in% colnames(spe)
  ) %>%
  group_by(donor) %>%
  summarise(
    n_xenium_cells = n(),
    n_valid_keys = sum(valid_key),
    n_matched_imc = sum(matched_imc),
    matched_imc_percent = round(100 * n_matched_imc / n_xenium_cells, 1),
    .groups = "drop"
  ) %>%
  arrange(factor(
    donor,
    levels = c("Validation1", "Validation2", "ExcludedValidation", "Validation3", "Validation4", "Validation5", "Validation6")
  ))

write.csv(donor_audit, file.path(output_dir, "seven_donor_join_audit.csv"), row.names = FALSE)

assay_audit <- data.frame(
  assay = assayNames(sfe7),
  matrix_class = vapply(
    assayNames(sfe7),
    function(nm) paste(class(assay(sfe7, nm)), collapse = "/"),
    character(1)
  ),
  n_features = nrow(sfe7),
  n_cells = ncol(sfe7),
  nonzero = vapply(
    assayNames(sfe7),
    function(nm) Matrix::nnzero(assay(sfe7, nm)),
    numeric(1)
  ),
  stringsAsFactors = FALSE
)
write.csv(assay_audit, file.path(output_dir, "seven_donor_assay_audit.csv"), row.names = FALSE)

message("Saving joined seven-donor SFE: ", joined_path)
saveRDS(sfe7, joined_path)

cat("\nJoined object\n")
show(sfe7)
cat("\nDonor audit\n")
print(donor_audit)
cat("\nAssay audit\n")
print(assay_audit)
cat("\nSaved:", normalizePath(joined_path, winslash = "/"), "\n")
