
## Transfer annotations:
##   (A) SFE -> SPE  (lv2/lv3, etc.)
##   (B) SPE -> SFE  (IMC clusters, etc.)
## ============================================
## ----------------------------
## USER INPUTS (edit these)
## ----------------------------
# Objects assumed already in memory:
#   spe  : SpatialExperiment (IMC)
#   sfe  : SpatialFeatureExperiment (Xenium / SFE)
#
# Required:
#   colData(sfe)$KEY exists and matches colnames(spe)

# Columns to transfer from SFE to SPE
SFE_TO_SPE_COLS <- c("lv2_anno", "lv3_anno")   # add more if you want

# Cluster column to transfer from SPE to SFE
SPE_CLUSTER_COL <- "Rphenograph_500"          # change if your column differs

# Output column names to create
SFE_TO_SPE_SUFFIX <- "_from_sfe"
SPE_TO_SFE_PREFIX <- "IMC_"


## ----------------------------
## Helper: check KEY overlap
## ----------------------------
check_key_overlap <- function(spe, sfe, key_col = "KEY") {
  spe_keys <- colnames(spe)
  sfe_keys <- as.character(colData(sfe)[[key_col]])
  n_overlap <- length(intersect(spe_keys, sfe_keys))
  message("KEY overlap (colnames(spe) vs colData(sfe)$", key_col, "): ", n_overlap)
  invisible(n_overlap)
}

## ----------------------------
## A) Transfer SFE -> SPE
## ----------------------------
transfer_sfe_to_spe <- function(spe, sfe, key_col = "KEY", cols = SFE_TO_SPE_COLS,
                                suffix = SFE_TO_SPE_SUFFIX) {
  
  if (!key_col %in% colnames(colData(sfe))) {
    stop("colData(sfe)$", key_col, " not found.")
  }
  
  sfe_map <- as.data.frame(colData(sfe))[, c(key_col, cols), drop = FALSE]
  sfe_map[[key_col]] <- as.character(sfe_map[[key_col]])
  sfe_map <- sfe_map[!duplicated(sfe_map[[key_col]]), , drop = FALSE]
  
  idx <- match(colnames(spe), sfe_map[[key_col]])
  
  for (cc in cols) {
    new_name <- paste0(cc, suffix)
    colData(spe)[[new_name]] <- sfe_map[[cc]][idx]
    message("Added to spe: colData(spe)$", new_name,
            " | non-NA = ", sum(!is.na(colData(spe)[[new_name]])))
  }
  
  spe
}

## ----------------------------
## B) Transfer SPE -> SFE
## ----------------------------
transfer_spe_to_sfe <- function(sfe, spe, key_col = "KEY", spe_cluster_col = SPE_CLUSTER_COL,
                                prefix = SPE_TO_SFE_PREFIX) {
  
  if (!key_col %in% colnames(colData(sfe))) {
    stop("colData(sfe)$", key_col, " not found.")
  }
  if (!spe_cluster_col %in% colnames(colData(spe))) {
    stop("colData(spe)$", spe_cluster_col, " not found.")
  }
  
  # Map cluster vector by SPE cell IDs (names = colnames(spe))
  imc_clust <- colData(spe)[[spe_cluster_col]]
  names(imc_clust) <- colnames(spe)
  
  idx <- match(as.character(colData(sfe)[[key_col]]), names(imc_clust))
  
  out_name <- paste0(prefix, spe_cluster_col)
  colData(sfe)[[out_name]] <- imc_clust[idx]
  
  message("Added to sfe: colData(sfe)$", out_name,
          " | non-NA = ", sum(!is.na(colData(sfe)[[out_name]])))
  
  sfe
}


## ----------------------------
## RUN (uncomment when ready)
## ----------------------------
# check_key_overlap(spe, sfe, key_col = "KEY")
# spe <- transfer_sfe_to_spe(spe, sfe, key_col = "KEY", cols = SFE_TO_SPE_COLS)
# sfe <- transfer_spe_to_sfe(sfe, spe, key_col = "KEY", spe_cluster_col = SPE_CLUSTER_COL)

## ----------------------------
## SANITY CHECKS (examples)
## ----------------------------
# table(is.na(colData(spe)$lv3_anno_from_sfe))
# table(is.na(colData(sfe)$IMC_Rphenograph_500), useNA = "ifany")
