source("data_preprocessing/common/config.R")
suppressPackageStartupMessages({
  library(SpatialFeatureExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(Seurat)
  library(scMetabolism)
  library(decoupleR)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

set.seed(20260729)
options(future.globals.maxSize = 8 * 1024^3)

project_root <- normalizePath(
  file.path(getwd(), "seven_donor_reanalysis"),
  winslash = "/",
  mustWork = TRUE
)
output_dir <- file.path(project_root, "outputs")
score_dir <- file.path(output_dir, "seven_donor_scores")
dir.create(score_dir, recursive = TRUE, showWarnings = FALSE)

sfe_path <- file.path(output_dir, "sfe_seven_donors_annotated_final.rds")
matrisome_path <- paste0(
  paste0(project_path("Revision/"), "/"),
  "Nan annotation/Data/Hs_Matrisome_Masterlist_Naba.csv"
)

scmetabolism_cache <- file.path(
  score_dir,
  "seven_donor_fib_ec_scmetabolism_score_matrix.rds"
)
ecm_cache <- file.path(
  score_dir,
  "seven_donor_fib_ec_matrisome_aucell.rds"
)

to_dgC <- function(x) {
  if (inherits(x, "dgCMatrix")) return(x)
  as(as(x, "CsparseMatrix"), "dgCMatrix")
}

message("Reading finalized seven-donor SFE")
sfe7 <- readRDS(sfe_path)
metadata7 <- as.data.frame(colData(sfe7), check.names = FALSE)
rownames(metadata7) <- colnames(sfe7)

keep <- as.character(metadata7$lv1_anno_final) %in% c("fibroblast", "endothelial")
sfe_stromal <- sfe7[, keep]
metadata_stromal <- as.data.frame(colData(sfe_stromal), check.names = FALSE)
rownames(metadata_stromal) <- colnames(sfe_stromal)

cell_manifest <- metadata_stromal %>%
  rownames_to_column("cell") %>%
  transmute(
    cell,
    donor = as.character(SampleId),
    lv1_anno = as.character(lv1_anno_final),
    lv2_anno = as.character(lv2_anno_final),
    metfiblabel = as.character(metfiblabel),
    metEClabel = as.character(metEClabel),
    KEY = as.character(KEY),
    imc_matched = as.logical(imc_matched)
  )
write.csv(
  cell_manifest,
  file.path(score_dir, "seven_donor_fib_ec_cell_manifest.csv"),
  row.names = FALSE
)

message("Running or reading seven-donor scMetabolism AUCell scores")
if (file.exists(scmetabolism_cache)) {
  scmetabolism_scores <- readRDS(scmetabolism_cache)
} else {
  counts <- to_dgC(assay(sfe_stromal, "counts"))

  old_assay_version <- getOption("Seurat.object.assay.version")
  options(Seurat.object.assay.version = "v3")
  seu <- CreateSeuratObject(
    counts = counts,
    meta.data = metadata_stromal,
    assay = "RNA",
    project = "seven_donor_fib_ec_scMetabolism"
  )
  if (is.null(old_assay_version)) {
    options(Seurat.object.assay.version = NULL)
  } else {
    options(Seurat.object.assay.version = old_assay_version)
  }

  seu <- NormalizeData(seu, verbose = FALSE)
  physical_cores <- parallel::detectCores(logical = FALSE)
  ncores <- if (is.na(physical_cores)) 1L else min(2L, physical_cores)

  seu <- scMetabolism::sc.metabolism.Seurat(
    obj = seu,
    method = "AUCell",
    imputation = FALSE,
    ncores = ncores,
    metabolism.type = "KEGG"
  )
  scmetabolism_scores <- seu@assays$METABOLISM$score
  saveRDS(scmetabolism_scores, scmetabolism_cache)
  rm(seu)
  gc()
}

required_pathways <- c(
  "Glycolysis / Gluconeogenesis",
  "Citrate cycle (TCA cycle)"
)
missing_pathways <- setdiff(required_pathways, rownames(scmetabolism_scores))
if (length(missing_pathways)) {
  stop(
    "Required scMetabolism pathways are missing: ",
    paste(missing_pathways, collapse = ", ")
  )
}

score_cells <- colnames(scmetabolism_scores)
if (!all(score_cells %in% cell_manifest$cell)) {
  seurat_to_original <- setNames(
    cell_manifest$cell,
    make.names(cell_manifest$cell, unique = TRUE)
  )
  restored_cells <- unname(seurat_to_original[score_cells])
  if (anyNA(restored_cells) || anyDuplicated(restored_cells)) {
    stop("Could not restore scMetabolism cell names to the joined SFE identifiers.")
  }
  colnames(scmetabolism_scores) <- restored_cells
  score_cells <- restored_cells
  saveRDS(scmetabolism_scores, scmetabolism_cache)
}
stopifnot(all(score_cells %in% cell_manifest$cell))
energy_score_wide <- as.data.frame(
  t(scmetabolism_scores[required_pathways, , drop = FALSE]),
  check.names = FALSE
)
energy_score_wide$oxphos_scmetabolism <- if (
    "Oxidative phosphorylation" %in% rownames(scmetabolism_scores)) {
  as.numeric(scmetabolism_scores[
    "Oxidative phosphorylation",
    rownames(energy_score_wide)
  ])
} else {
  NA_real_
}

energy_scores <- energy_score_wide %>%
  rownames_to_column("cell") %>%
  left_join(cell_manifest, by = "cell") %>%
  transmute(
    cell,
    donor,
    lv1_anno,
    lv2_anno,
    metfiblabel,
    metEClabel,
    KEY,
    imc_matched,
    glycolysis_scmetabolism = .data[["Glycolysis / Gluconeogenesis"]],
    tca_scmetabolism = .data[["Citrate cycle (TCA cycle)"]],
    oxphos_scmetabolism
  )

write.csv(
  energy_scores,
  file.path(score_dir, "seven_donor_scmetabolism_energy_cell_scores.csv"),
  row.names = FALSE
)

message("Running or reading seven-donor matrisome AUCell scores")
matrisome <- read.csv(
  matrisome_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
ecm_core <- matrisome %>%
  filter(`Matrisome Division` == "Core matrisome") %>%
  transmute(
    source = as.character(`Matrisome Category`),
    target = as.character(`Gene Symbol`)
  )
ecm_network <- bind_rows(
  ecm_core,
  ecm_core %>% transmute(source = "CoreMatrisome", target)
) %>%
  distinct(source, target)

gene_availability <- ecm_network %>%
  group_by(source) %>%
  summarise(
    n_genes = n_distinct(target),
    n_detected = sum(unique(target) %in% rownames(sfe_stromal)),
    detected_fraction = n_detected / n_genes,
    .groups = "drop"
  )
write.csv(
  gene_availability,
  file.path(score_dir, "seven_donor_matrisome_gene_availability.csv"),
  row.names = FALSE
)

if (file.exists(ecm_cache)) {
  ecm_scores <- readRDS(ecm_cache)
} else {
  ecm_scores <- decoupleR::run_aucell(
    mat = to_dgC(assay(sfe_stromal, "logcounts")),
    network = ecm_network,
    minsize = 0,
    .source = source,
    .target = target,
    nproc = 1
  )
  saveRDS(ecm_scores, ecm_cache)
}

ecm_score_long <- ecm_scores %>%
  transmute(
    ecm_score = paste0(as.character(source), " score"),
    cell = as.character(condition),
    score = as.numeric(score)
  ) %>%
  left_join(cell_manifest, by = "cell")

write.csv(
  ecm_score_long,
  file.path(score_dir, "seven_donor_matrisome_aucell_cell_scores.csv"),
  row.names = FALSE
)

audit <- tibble(
  metric = c(
    "Joined SFE cells",
    "Annotated fibroblast cells scored",
    "Annotated endothelial cells scored",
    "scMetabolism pathways",
    "Matrisome score sets"
  ),
  value = c(
    ncol(sfe7),
    sum(cell_manifest$lv1_anno == "fibroblast"),
    sum(cell_manifest$lv1_anno == "endothelial"),
    nrow(scmetabolism_scores),
    n_distinct(ecm_score_long$ecm_score)
  )
)
write.csv(
  audit,
  file.path(score_dir, "seven_donor_score_audit.csv"),
  row.names = FALSE
)

cat("\nSeven-donor score audit\n")
print(audit)
cat("\nSaved scores in:", normalizePath(score_dir, winslash = "/"), "\n")
