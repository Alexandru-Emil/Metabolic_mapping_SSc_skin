source("data_preprocessing/common/config.R")
suppressPackageStartupMessages({
  library(SpatialFeatureExperiment)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(Matrix)
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(RANN)
  library(mclust)
  library(aricode)
})

set.seed(20260729)
options(future.globals.maxSize = 8 * 1024^3)

project_root <- normalizePath(
  file.path(getwd(), "seven_donor_reanalysis"),
  winslash = "/",
  mustWork = TRUE
)
output_dir <- file.path(project_root, "outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

joined_path <- file.path(output_dir, "sfe_seven_donors_joined_unannotated.rds")
annotated_path <- file.path(output_dir, "sfe_seven_donors_annotated_final.rds")

reference_seurat_path <- paste0(
  paste0(project_path("Revision/"), "/"),
  "Nan annotation/Data/Xenium sfe objects/",
  "VedaSScSkinXe_EcFibAnno_260225234740.rds"
)
spe_path <- paste0(
  paste0(project_path("Revision/"), "/"),
  "Nan annotation/Data/IMC spe objects/New IMC/spe2025-12-17 umap.RData"
)

donor_levels <- c("Validation1", "Validation2", "ExcludedValidation", "Validation3", "Validation4", "Validation5", "Validation6")
major_levels <- c(
  "epithelial",
  "fibroblast",
  "pericyte",
  "endothelial",
  "muscle",
  "telocyte",
  "Schwann cell",
  "myeloid cell",
  "lymphocyte"
)
fib_levels <- c(
  "papillary_Fib",
  "COMP_Fib",
  "COL8A1_Fib",
  "PI16_Fib",
  "CXCL12_Fib",
  "CCL19_Fib",
  "NGFR_Fib",
  "COCH_Fib"
)
ec_levels <- c(
  "cycling_EC",
  "ACKR1_EC",
  "EPC",
  "ACTA2_EC",
  "HEY1_EC",
  "lymphatic_EC"
)

to_dgC <- function(x) {
  if (inherits(x, "dgCMatrix")) return(x)
  as(as(x, "CsparseMatrix"), "dgCMatrix")
}

safe_ari <- function(x, y) {
  keep <- !is.na(x) & !is.na(y)
  if (sum(keep) < 2L) return(NA_real_)
  mclust::adjustedRandIndex(as.character(x[keep]), as.character(y[keep]))
}

safe_nmi <- function(x, y) {
  keep <- !is.na(x) & !is.na(y)
  if (sum(keep) < 2L) return(NA_real_)
  aricode::NMI(as.character(x[keep]), as.character(y[keep]))
}

majority_cluster_map <- function(clusters, labels, reference_mask, donors) {
  clusters <- as.character(clusters)
  labels <- as.character(labels)
  donors <- as.character(donors)
  reference_mask <- reference_mask & !is.na(labels) & nzchar(labels)

  tab <- table(clusters[reference_mask], labels[reference_mask])
  mapping <- if (nrow(tab)) {
    apply(tab, 1, function(x) names(which.max(x)))
  } else {
    character()
  }
  predictions <- unname(mapping[clusters])

  cluster_names <- sort(unique(clusters))
  cluster_stats <- bind_rows(lapply(cluster_names, function(cluster_id) {
    all_idx <- clusters == cluster_id
    ref_idx <- all_idx & reference_mask
    label_tab <- sort(table(labels[ref_idx]), decreasing = TRUE)
    n_ref <- sum(ref_idx)
    top_label <- if (length(label_tab)) names(label_tab)[1] else NA_character_
    top_n <- if (length(label_tab)) as.integer(label_tab[1]) else 0L
    data.frame(
      cluster = cluster_id,
      assigned_label = top_label,
      n_all = sum(all_idx),
      n_reference = n_ref,
      n_ssc230 = sum(all_idx & donors == "ExcludedValidation"),
      n_donors = length(unique(donors[all_idx])),
      n_reference_donors = length(unique(donors[ref_idx])),
      purity = if (n_ref > 0L) top_n / n_ref else NA_real_,
      stringsAsFactors = FALSE
    )
  }))

  list(predictions = predictions, mapping = mapping, cluster_stats = cluster_stats)
}

evaluate_cluster_column <- function(
    metadata,
    cluster_col,
    reference_label_col,
    query_donor = "ExcludedValidation",
    target_clusters) {
  clusters <- as.character(metadata[[cluster_col]])
  labels <- as.character(metadata[[reference_label_col]])
  donors <- as.character(metadata[["SampleId"]])
  reference_mask <- donors != query_donor & !is.na(labels) & nzchar(labels)

  mapped <- majority_cluster_map(clusters, labels, reference_mask, donors)
  pred <- mapped$predictions
  query_clusters <- unique(clusters[donors == query_donor])
  query_stats <- mapped$cluster_stats %>%
    filter(cluster %in% query_clusters)

  accuracy <- mean(pred[reference_mask] == labels[reference_mask], na.rm = TRUE)
  mean_donors <- if (nrow(query_stats)) mean(query_stats$n_donors) else 0
  min_ref_support <- if (nrow(query_stats)) min(query_stats$n_reference) else 0
  unsupported <- if (nrow(query_stats)) sum(query_stats$n_reference == 0) else 0
  n_clusters <- length(unique(clusters))
  ari <- safe_ari(clusters[reference_mask], labels[reference_mask])
  nmi <- safe_nmi(clusters[reference_mask], labels[reference_mask])

  selection_score <- accuracy +
    0.05 * ifelse(is.finite(ari), ari, 0) +
    0.03 * pmin(mean_donors / 7, 1) -
    0.003 * abs(n_clusters - target_clusters) -
    0.05 * unsupported

  data.frame(
    cluster_col = cluster_col,
    n_clusters = n_clusters,
    reference_label_accuracy = accuracy,
    adjusted_rand_index = ari,
    normalized_mutual_information = nmi,
    mean_donors_in_ssc230_clusters = mean_donors,
    min_reference_cells_in_ssc230_clusters = min_ref_support,
    unsupported_ssc230_clusters = unsupported,
    selection_score = selection_score,
    stringsAsFactors = FALSE
  )
}

add_resolution_grid <- function(
    object,
    graph_name,
    resolutions,
    prefix,
    reference_label_col,
    target_clusters) {
  metrics <- vector("list", length(resolutions))

  for (i in seq_along(resolutions)) {
    resolution <- resolutions[[i]]
    cluster_col <- paste0(prefix, "_res_", gsub("\\.", "_", format(resolution, trim = TRUE)))
    object <- FindClusters(
      object,
      graph.name = graph_name,
      resolution = resolution,
      algorithm = 4,
      random.seed = 20260729,
      verbose = FALSE
    )
    object[[cluster_col]] <- as.character(Idents(object))
    metrics[[i]] <- evaluate_cluster_column(
      object@meta.data,
      cluster_col = cluster_col,
      reference_label_col = reference_label_col,
      target_clusters = target_clusters
    ) %>%
      mutate(resolution = resolution, .before = 1)
  }

  metrics <- bind_rows(metrics) %>%
    mutate(
      eligible = unsupported_ssc230_clusters == 0 &
        min_reference_cells_in_ssc230_clusters >= 20 &
        mean_donors_in_ssc230_clusters >= 4
    )

  candidate <- metrics %>% filter(eligible)
  if (!nrow(candidate)) candidate <- metrics
  selected <- candidate %>%
    arrange(desc(selection_score), abs(n_clusters - target_clusters), resolution) %>%
    slice(1)

  list(
    object = object,
    metrics = metrics,
    selected_col = selected$cluster_col[[1]],
    selected_resolution = selected$resolution[[1]]
  )
}

knn_label <- function(reference_embedding, reference_labels, query_embedding, k = 25L) {
  reference_labels <- as.character(reference_labels)
  keep <- !is.na(reference_labels) & nzchar(reference_labels)
  reference_embedding <- reference_embedding[keep, , drop = FALSE]
  reference_labels <- reference_labels[keep]

  if (!nrow(query_embedding)) return(character())
  k <- min(k, nrow(reference_embedding))
  neighbors <- RANN::nn2(reference_embedding, query_embedding, k = k)$nn.idx
  apply(neighbors, 1, function(idx) {
    tab <- sort(table(reference_labels[idx]), decreasing = TRUE)
    names(tab)[1]
  })
}

run_harmony_model <- function(
    counts,
    metadata,
    project,
    nfeatures = 3000L,
    npcs = 30L,
    umap_name = "umap_harmony") {
  metadata <- as.data.frame(metadata, check.names = FALSE)
  rownames(metadata) <- colnames(counts)
  metadata$SampleId <- factor(as.character(metadata$SampleId), levels = donor_levels)

  old_assay_version <- getOption("Seurat.object.assay.version")
  options(Seurat.object.assay.version = "v3")
  object <- CreateSeuratObject(
    counts = to_dgC(counts),
    meta.data = metadata,
    assay = "RNA",
    project = project
  )
  if (is.null(old_assay_version)) {
    options(Seurat.object.assay.version = NULL)
  } else {
    options(Seurat.object.assay.version = old_assay_version)
  }

  object <- NormalizeData(
    object,
    normalization.method = "LogNormalize",
    scale.factor = 1000,
    verbose = FALSE
  )
  object <- FindVariableFeatures(
    object,
    selection.method = "vst",
    nfeatures = min(nfeatures, nrow(object)),
    verbose = FALSE
  )
  object <- ScaleData(object, features = VariableFeatures(object), verbose = FALSE)
  object <- RunPCA(
    object,
    features = VariableFeatures(object),
    npcs = npcs,
    seed.use = 20260729,
    verbose = FALSE
  )
  object <- harmony::RunHarmony(
    object,
    group.by.vars = "SampleId",
    reduction.use = "pca",
    dims.use = seq_len(npcs),
    reduction.save = "harmony",
    project.dim = FALSE,
    plot_convergence = FALSE,
    verbose = FALSE
  )
  object <- RunUMAP(
    object,
    reduction = "harmony",
    dims = seq_len(npcs),
    reduction.name = umap_name,
    reduction.key = paste0(gsub("[^A-Za-z0-9]", "", umap_name), "_"),
    seed.use = 20260729,
    verbose = FALSE
  )
  object <- FindNeighbors(
    object,
    reduction = "harmony",
    dims = seq_len(npcs),
    graph.name = c("harmony_nn", "harmony_snn"),
    verbose = FALSE
  )
  object
}

extract_joint_labels <- function(object, cluster_col, reference_label_col) {
  metadata <- object@meta.data
  clusters <- as.character(metadata[[cluster_col]])
  ref_labels <- as.character(metadata[[reference_label_col]])
  donors <- as.character(metadata$SampleId)
  reference_mask <- donors != "ExcludedValidation" & !is.na(ref_labels) & nzchar(ref_labels)

  mapped <- majority_cluster_map(clusters, ref_labels, reference_mask, donors)
  predictions <- mapped$predictions

  missing_query <- donors == "ExcludedValidation" & (is.na(predictions) | !nzchar(predictions))
  if (any(missing_query)) {
    embedding <- Embeddings(object, "harmony")
    predictions[missing_query] <- knn_label(
      embedding[reference_mask, , drop = FALSE],
      ref_labels[reference_mask],
      embedding[missing_query, , drop = FALSE],
      k = 25L
    )
  }

  list(labels = predictions, cluster_stats = mapped$cluster_stats)
}

run_joint_subcluster <- function(
    counts_all,
    metadata_all,
    parent_joint_labels,
    parent_label,
    reference_subtype_labels,
    subtype_levels,
    target_clusters,
  file_prefix) {
  donors <- as.character(metadata_all$SampleId)
  reference_parent <- donors != "ExcludedValidation" &
    !is.na(reference_subtype_labels) &
    reference_subtype_labels %in% subtype_levels
  query_parent <- donors == "ExcludedValidation" & parent_joint_labels == parent_label
  keep <- reference_parent | query_parent

  counts <- counts_all[, keep, drop = FALSE]
  metadata <- metadata_all[keep, , drop = FALSE]
  metadata$reference_subtype <- as.character(reference_subtype_labels[keep])

  object <- run_harmony_model(
    counts,
    metadata,
    project = paste0("seven_donor_", file_prefix),
    nfeatures = 2500L,
    npcs = 30L,
    umap_name = paste0("umap_", file_prefix, "_harmony")
  )
  grid <- add_resolution_grid(
    object,
    graph_name = "harmony_snn",
    resolutions = c(0.2, 0.4, 0.6, 0.8, 1, 1.2),
    prefix = paste0("joint_", file_prefix),
    reference_label_col = "reference_subtype",
    target_clusters = target_clusters
  )
  object <- grid$object
  selected <- extract_joint_labels(
    object,
    cluster_col = grid$selected_col,
    reference_label_col = "reference_subtype"
  )

  subtype <- selected$labels
  subtype[!subtype %in% subtype_levels] <- NA_character_
  embedding <- Embeddings(object, paste0("umap_", file_prefix, "_harmony"))

  result <- data.frame(
    cell = colnames(object),
    SampleId = as.character(object$SampleId),
    reference_subtype = as.character(object$reference_subtype),
    joint_subtype = subtype,
    joint_cluster = as.character(object@meta.data[[grid$selected_col]]),
    UMAP1 = embedding[, 1],
    UMAP2 = embedding[, 2],
    stringsAsFactors = FALSE
  )

  write.csv(
    grid$metrics,
    file.path(output_dir, paste0(file_prefix, "_joint_resolution_grid.csv")),
    row.names = FALSE
  )
  write.csv(
    selected$cluster_stats,
    file.path(output_dir, paste0(file_prefix, "_joint_cluster_label_map.csv")),
    row.names = FALSE
  )
  write.csv(
    result,
    file.path(output_dir, paste0(file_prefix, "_joint_labels_and_umap.csv")),
    row.names = FALSE
  )

  rm(object)
  gc()

  list(
    result = result,
    metrics = grid$metrics,
    selected_col = grid$selected_col,
    selected_resolution = grid$selected_resolution
  )
}

run_anchor_transfer <- function(
    reference_counts,
    reference_labels,
    query_counts,
    nfeatures = 3000L,
    npcs = 30L,
    project = "anchor_transfer") {
  reference_labels <- as.character(reference_labels)
  keep_reference <- !is.na(reference_labels) & nzchar(reference_labels)
  reference_counts <- reference_counts[, keep_reference, drop = FALSE]
  reference_labels <- reference_labels[keep_reference]

  if (!ncol(query_counts)) {
    return(data.frame())
  }

  ref_meta <- data.frame(
    transfer_label = reference_labels,
    row.names = colnames(reference_counts),
    stringsAsFactors = FALSE
  )
  query_meta <- data.frame(
    row.names = colnames(query_counts),
    stringsAsFactors = FALSE
  )

  old_assay_version <- getOption("Seurat.object.assay.version")
  options(Seurat.object.assay.version = "v3")
  reference <- CreateSeuratObject(
    counts = to_dgC(reference_counts),
    meta.data = ref_meta,
    assay = "RNA",
    project = paste0(project, "_reference")
  )
  query <- CreateSeuratObject(
    counts = to_dgC(query_counts),
    meta.data = query_meta,
    assay = "RNA",
    project = paste0(project, "_query")
  )
  if (is.null(old_assay_version)) {
    options(Seurat.object.assay.version = NULL)
  } else {
    options(Seurat.object.assay.version = old_assay_version)
  }

  reference <- NormalizeData(reference, scale.factor = 1000, verbose = FALSE)
  query <- NormalizeData(query, scale.factor = 1000, verbose = FALSE)
  reference <- FindVariableFeatures(
    reference,
    selection.method = "vst",
    nfeatures = min(nfeatures, nrow(reference)),
    verbose = FALSE
  )
  reference <- ScaleData(
    reference,
    features = VariableFeatures(reference),
    verbose = FALSE
  )
  reference <- RunPCA(
    reference,
    features = VariableFeatures(reference),
    npcs = npcs,
    seed.use = 20260729,
    verbose = FALSE
  )

  k_filter <- if (ncol(query) > 25L) min(100L, ncol(query) - 1L) else NA
  anchors <- FindTransferAnchors(
    reference = reference,
    query = query,
    normalization.method = "LogNormalize",
    reference.reduction = "pca",
    reduction = "pcaproject",
    features = VariableFeatures(reference),
    dims = seq_len(npcs),
    k.filter = k_filter,
    verbose = FALSE
  )
  predictions <- TransferData(
    anchorset = anchors,
    refdata = reference$transfer_label,
    dims = seq_len(npcs),
    verbose = FALSE
  )
  predictions$cell <- rownames(predictions)
  predictions <- predictions %>% relocate(cell)

  rm(reference, query, anchors)
  gc()
  predictions
}

message("Reading joined seven-donor SFE")
sfe7 <- readRDS(joined_path)
counts7 <- to_dgC(assay(sfe7, "counts"))
metadata7 <- as.data.frame(colData(sfe7), check.names = FALSE)
rownames(metadata7) <- colnames(sfe7)
metadata7$SampleId <- as.character(metadata7$SampleId)

message("Reading finalized six-donor ResolVI annotation reference")
reference_seurat <- readRDS(reference_seurat_path)
reference_cell_key <- paste(
  as.character(reference_seurat$SampleId),
  colnames(reference_seurat),
  sep = "_"
)
stopifnot(!anyDuplicated(reference_cell_key))
reference_lv1 <- setNames(as.character(reference_seurat$lv1_anno), reference_cell_key)
reference_lv2 <- setNames(as.character(reference_seurat$lv2_anno), reference_cell_key)
rm(reference_seurat)
gc()

metadata7$reference_lv1 <- unname(reference_lv1[rownames(metadata7)])
metadata7$reference_lv2 <- unname(reference_lv2[rownames(metadata7)])

stopifnot(
  all(!is.na(metadata7$reference_lv1[metadata7$SampleId != "ExcludedValidation"])),
  all(is.na(metadata7$reference_lv1[metadata7$SampleId == "ExcludedValidation"]))
)

message("Running seven-donor joint Harmony integration")
joint_full <- run_harmony_model(
  counts7,
  metadata7,
  project = "seven_donor_joint_major",
  nfeatures = 3000L,
  npcs = 30L,
  umap_name = "umap_joint_harmony"
)

joint_grid <- add_resolution_grid(
  joint_full,
  graph_name = "harmony_snn",
  resolutions = c(0.2, 0.4, 0.6, 0.8, 1, 1.2),
  prefix = "joint_major",
  reference_label_col = "reference_lv1",
  target_clusters = 20L
)
joint_full <- joint_grid$object
joint_major <- extract_joint_labels(
  joint_full,
  cluster_col = joint_grid$selected_col,
  reference_label_col = "reference_lv1"
)
joint_major_labels <- joint_major$labels
joint_major_labels[!joint_major_labels %in% major_levels] <- NA_character_

joint_umap <- Embeddings(joint_full, "umap_joint_harmony")
joint_full_result <- data.frame(
  cell = colnames(joint_full),
  SampleId = as.character(joint_full$SampleId),
  reference_lv1 = as.character(joint_full$reference_lv1),
  joint_lv1 = joint_major_labels,
  joint_cluster = as.character(joint_full@meta.data[[joint_grid$selected_col]]),
  UMAP1 = joint_umap[, 1],
  UMAP2 = joint_umap[, 2],
  stringsAsFactors = FALSE
)

write.csv(
  joint_grid$metrics,
  file.path(output_dir, "major_joint_resolution_grid.csv"),
  row.names = FALSE
)
write.csv(
  joint_major$cluster_stats,
  file.path(output_dir, "major_joint_cluster_label_map.csv"),
  row.names = FALSE
)
write.csv(
  joint_full_result,
  file.path(output_dir, "major_joint_labels_and_umap.csv"),
  row.names = FALSE
)

rm(joint_full)
gc()

message("Running seven-donor joint fibroblast re-clustering")
joint_fib <- run_joint_subcluster(
  counts_all = counts7,
  metadata_all = metadata7,
  parent_joint_labels = joint_major_labels,
  parent_label = "fibroblast",
  reference_subtype_labels = metadata7$reference_lv2,
  subtype_levels = fib_levels,
  target_clusters = 9L,
  file_prefix = "fibroblast"
)

message("Running seven-donor joint EC re-clustering")
joint_ec <- run_joint_subcluster(
  counts_all = counts7,
  metadata_all = metadata7,
  parent_joint_labels = joint_major_labels,
  parent_label = "endothelial",
  reference_subtype_labels = metadata7$reference_lv2,
  subtype_levels = ec_levels,
  target_clusters = 7L,
  file_prefix = "endothelial"
)

joint_lv2 <- setNames(rep(NA_character_, ncol(sfe7)), colnames(sfe7))
joint_lv2[joint_fib$result$cell] <- joint_fib$result$joint_subtype
joint_lv2[joint_ec$result$cell] <- joint_ec$result$joint_subtype

message("Running independent six-donor reference mapping: major populations")
reference_cells <- metadata7$SampleId != "ExcludedValidation"
query_cells <- metadata7$SampleId == "ExcludedValidation"

map_major <- run_anchor_transfer(
  reference_counts = counts7[, reference_cells, drop = FALSE],
  reference_labels = metadata7$reference_lv1[reference_cells],
  query_counts = counts7[, query_cells, drop = FALSE],
  nfeatures = 3000L,
  npcs = 30L,
  project = "ssc230_major"
)
colnames(map_major)[colnames(map_major) == "predicted.id"] <- "mapped_lv1"
colnames(map_major)[colnames(map_major) == "prediction.score.max"] <- "mapped_lv1_score"
write.csv(map_major, file.path(output_dir, "ssc230_major_reference_mapping.csv"), row.names = FALSE)

mapped_lv1_query <- setNames(map_major$mapped_lv1, map_major$cell)

message("Running independent reference mapping: fibroblast subtypes")
query_fib_cells <- names(mapped_lv1_query)[mapped_lv1_query == "fibroblast"]
reference_fib_cells <- rownames(metadata7)[
  reference_cells & metadata7$reference_lv1 == "fibroblast" &
    !is.na(metadata7$reference_lv2)
]
map_fib <- run_anchor_transfer(
  reference_counts = counts7[, reference_fib_cells, drop = FALSE],
  reference_labels = metadata7$reference_lv2[match(reference_fib_cells, rownames(metadata7))],
  query_counts = counts7[, query_fib_cells, drop = FALSE],
  nfeatures = 2500L,
  npcs = 30L,
  project = "ssc230_fibroblast"
)
if (nrow(map_fib)) {
  colnames(map_fib)[colnames(map_fib) == "predicted.id"] <- "mapped_lv2"
  colnames(map_fib)[colnames(map_fib) == "prediction.score.max"] <- "mapped_lv2_score"
}
write.csv(map_fib, file.path(output_dir, "ssc230_fibroblast_reference_mapping.csv"), row.names = FALSE)

message("Running independent reference mapping: EC subtypes")
query_ec_cells <- names(mapped_lv1_query)[mapped_lv1_query == "endothelial"]
reference_ec_cells <- rownames(metadata7)[
  reference_cells & metadata7$reference_lv1 == "endothelial" &
    !is.na(metadata7$reference_lv2)
]
map_ec <- run_anchor_transfer(
  reference_counts = counts7[, reference_ec_cells, drop = FALSE],
  reference_labels = metadata7$reference_lv2[match(reference_ec_cells, rownames(metadata7))],
  query_counts = counts7[, query_ec_cells, drop = FALSE],
  nfeatures = 2500L,
  npcs = 30L,
  project = "ssc230_endothelial"
)
if (nrow(map_ec)) {
  colnames(map_ec)[colnames(map_ec) == "predicted.id"] <- "mapped_lv2"
  colnames(map_ec)[colnames(map_ec) == "prediction.score.max"] <- "mapped_lv2_score"
}
write.csv(map_ec, file.path(output_dir, "ssc230_endothelial_reference_mapping.csv"), row.names = FALSE)

mapped_lv2_query <- c(
  if (nrow(map_fib)) setNames(map_fib$mapped_lv2, map_fib$cell) else character(),
  if (nrow(map_ec)) setNames(map_ec$mapped_lv2, map_ec$cell) else character()
)
mapped_lv2_score_query <- c(
  if (nrow(map_fib)) setNames(map_fib$mapped_lv2_score, map_fib$cell) else numeric(),
  if (nrow(map_ec)) setNames(map_ec$mapped_lv2_score, map_ec$cell) else numeric()
)

joint_lv1_named <- setNames(joint_major_labels, colnames(sfe7))
mapping_lv1_all <- metadata7$reference_lv1
mapping_lv1_all[query_cells] <- unname(mapped_lv1_query[rownames(metadata7)[query_cells]])
mapping_lv2_all <- metadata7$reference_lv2
mapping_lv2_all[query_cells] <- unname(mapped_lv2_query[rownames(metadata7)[query_cells]])

joint_lv2_named <- joint_lv2[colnames(sfe7)]

query_names <- rownames(metadata7)[query_cells]
comparison_query <- data.frame(
  cell = query_names,
  joint_lv1 = unname(joint_lv1_named[query_names]),
  mapped_lv1 = unname(mapping_lv1_all[query_cells]),
  mapped_lv1_score = unname(setNames(map_major$mapped_lv1_score, map_major$cell)[query_names]),
  joint_lv2 = unname(joint_lv2_named[query_names]),
  mapped_lv2 = unname(mapping_lv2_all[query_cells]),
  mapped_lv2_score = unname(mapped_lv2_score_query[query_names]),
  stringsAsFactors = FALSE
)
comparison_query$lv1_agree <- comparison_query$joint_lv1 == comparison_query$mapped_lv1
comparison_query$lv2_agree <- comparison_query$joint_lv2 == comparison_query$mapped_lv2

reference_names <- rownames(metadata7)[reference_cells]
joint_reference_lv1_accuracy <- mean(
  joint_lv1_named[reference_names] == metadata7$reference_lv1[reference_cells],
  na.rm = TRUE
)
joint_reference_lv2_keep <- !is.na(metadata7$reference_lv2[reference_cells])
joint_reference_lv2_accuracy <- mean(
  joint_lv2_named[reference_names][joint_reference_lv2_keep] ==
    metadata7$reference_lv2[reference_cells][joint_reference_lv2_keep],
  na.rm = TRUE
)
query_lv1_agreement <- mean(comparison_query$lv1_agree, na.rm = TRUE)
query_lv2_comparable <- !is.na(comparison_query$joint_lv2) &
  !is.na(comparison_query$mapped_lv2)
query_lv2_agreement <- mean(
  comparison_query$lv2_agree[query_lv2_comparable],
  na.rm = TRUE
)

comparison_summary <- data.frame(
  metric = c(
    "Six-donor major-label retention in joint clustering",
    "Six-donor Fib/EC subtype retention in joint clustering",
    "ExcludedValidation major-label agreement: joint versus mapping",
    "ExcludedValidation Fib/EC subtype agreement: joint versus mapping",
    "ExcludedValidation cells compared for major labels",
    "ExcludedValidation cells compared for Fib/EC subtypes"
  ),
  value = c(
    joint_reference_lv1_accuracy,
    joint_reference_lv2_accuracy,
    query_lv1_agreement,
    query_lv2_agreement,
    nrow(comparison_query),
    sum(query_lv2_comparable)
  ),
  unit = c("proportion", "proportion", "proportion", "proportion", "cells", "cells"),
  stringsAsFactors = FALSE
)

comparable <- isTRUE(
  joint_reference_lv1_accuracy >= 0.85 &&
    joint_reference_lv2_accuracy >= 0.75 &&
    query_lv1_agreement >= 0.85 &&
    query_lv2_agreement >= 0.70
)
recommended_strategy <- if (comparable) "reference_mapping" else "joint_harmony"

decision <- data.frame(
  recommended_strategy = recommended_strategy,
  comparable_by_predefined_thresholds = comparable,
  rationale = if (comparable) {
    paste(
      "The strategies met the predefined concordance thresholds. Reference mapping",
      "is preferred because it preserves every finalized six-donor ResolVI label",
      "and adds ExcludedValidation without retroactively redefining the published annotation space."
    )
  } else {
    paste(
      "The strategies did not meet all predefined concordance thresholds. The joint",
      "Harmony annotation is retained as the primary candidate pending marker-level review."
    )
  },
  stringsAsFactors = FALSE
)

write.csv(
  comparison_query,
  file.path(output_dir, "ssc230_joint_vs_mapping_cell_concordance.csv"),
  row.names = FALSE
)
write.csv(
  comparison_summary,
  file.path(output_dir, "annotation_strategy_comparison_summary.csv"),
  row.names = FALSE
)
write.csv(
  decision,
  file.path(output_dir, "annotation_strategy_decision.csv"),
  row.names = FALSE
)
write.csv(
  as.data.frame(table(
    joint = comparison_query$joint_lv1,
    mapped = comparison_query$mapped_lv1
  )),
  file.path(output_dir, "ssc230_major_label_crosstab.csv"),
  row.names = FALSE
)
write.csv(
  as.data.frame(table(
    joint = comparison_query$joint_lv2,
    mapped = comparison_query$mapped_lv2
  )),
  file.path(output_dir, "ssc230_subtype_label_crosstab.csv"),
  row.names = FALSE
)

message("Attaching both candidate annotations and selected final labels to the SFE")
sfe7$lv1_anno_joint <- unname(joint_lv1_named[colnames(sfe7)])
sfe7$lv2_anno_joint <- unname(joint_lv2_named[colnames(sfe7)])
sfe7$lv1_anno_mapped <- unname(mapping_lv1_all)
sfe7$lv2_anno_mapped <- unname(mapping_lv2_all)

mapped_major_score <- setNames(map_major$mapped_lv1_score, map_major$cell)
sfe7$lv1_mapping_score <- NA_real_
sfe7$lv1_mapping_score[query_cells] <- unname(mapped_major_score[query_names])
sfe7$lv2_mapping_score <- NA_real_
sfe7$lv2_mapping_score[query_cells] <- unname(mapped_lv2_score_query[query_names])

if (recommended_strategy == "reference_mapping") {
  sfe7$lv1_anno_final <- sfe7$lv1_anno_mapped
  sfe7$lv2_anno_final <- sfe7$lv2_anno_mapped
} else {
  sfe7$lv1_anno_final <- sfe7$lv1_anno_joint
  sfe7$lv2_anno_final <- sfe7$lv2_anno_joint
}

sfe7$lv1_anno_final <- factor(as.character(sfe7$lv1_anno_final), levels = major_levels)
sfe7$lv2_anno_final <- factor(
  as.character(sfe7$lv2_anno_final),
  levels = c(fib_levels, ec_levels)
)
sfe7$lv1_anno_new_text <- sfe7$lv1_anno_final
sfe7$lv2_anno_new <- sfe7$lv2_anno_final
sfe7$lv3_anno <- sfe7$lv2_anno_final
sfe7$annotation_strategy_final <- recommended_strategy

load(spe_path)
stopifnot(exists("spe"))
spe$metfiblabel <- case_when(
  as.character(spe$Rphenograph_500) %in% c("1", "2") ~ "Met_hi_Fib",
  as.character(spe$Rphenograph_500) == "3" ~ "Other_Fib",
  TRUE ~ NA_character_
)
spe$metEClabel <- case_when(
  as.character(spe$Rphenograph_200EC) %in% c("1", "2") ~ "Met_hi_EC",
  as.character(spe$Rphenograph_200EC) == "3" ~ "Other_EC",
  TRUE ~ NA_character_
)
key_idx <- match(as.character(sfe7$KEY), colnames(spe))
sfe7$metfiblabel <- as.character(spe$metfiblabel[key_idx])
sfe7$metEClabel <- as.character(spe$metEClabel[key_idx])
sfe7$imc_matched <- !is.na(key_idx)

metadata(sfe7)$annotation_comparison <- list(
  joint_method = paste(
    "Seven-donor LogNormalize/PCA/Harmony integration; Leiden clustering over",
    "a resolution grid; cluster labels anchored by majority finalized six-donor",
    "ResolVI annotations; fibroblasts and ECs re-integrated and re-clustered."
  ),
  mapping_method = paste(
    "Independent two-stage Seurat RNA anchor transfer from the finalized",
    "six-donor reference: major population, then fibroblast or EC subtype."
  ),
  selected_strategy = recommended_strategy,
  predefined_thresholds = c(
    joint_reference_major_retention = 0.85,
    joint_reference_subtype_retention = 0.75,
    ssc230_major_agreement = 0.85,
    ssc230_subtype_agreement = 0.70
  ),
  comparison_summary = comparison_summary,
  decision = decision
)

saveRDS(sfe7, annotated_path)

final_counts <- as.data.frame(table(
  donor = as.character(sfe7$SampleId),
  lv1_anno_final = as.character(sfe7$lv1_anno_final)
)) %>%
  filter(Freq > 0) %>%
  rename(n_cells = Freq)
write.csv(final_counts, file.path(output_dir, "seven_donor_final_major_label_counts.csv"), row.names = FALSE)

cat("\nAnnotation strategy comparison\n")
print(comparison_summary)
cat("\nDecision\n")
print(decision)
cat("\nSelected major resolution:", joint_grid$selected_resolution, "\n")
cat("Selected fibroblast resolution:", joint_fib$selected_resolution, "\n")
cat("Selected EC resolution:", joint_ec$selected_resolution, "\n")
cat("Saved annotated SFE:", normalizePath(annotated_path, winslash = "/"), "\n")
