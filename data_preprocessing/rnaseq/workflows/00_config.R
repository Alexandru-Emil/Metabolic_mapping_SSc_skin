# Scientific workflow derived from Revision/RNAseq PCS met inhibitors/results/scripts/00_config.R
# Run from the repository root with METABOLIC_INPUT_DIR set to the external project data.
source("data_preprocessing/common/config.R")
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(sva)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(openxlsx)
  library(rmarkdown)
  library(knitr)
})

set.seed(42)

pipeline_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

script_dir <- pipeline_script_dir()
root_dir <- normalizePath(
  Sys.getenv("RNASEQ_RESULTS_ROOT", dirname(script_dir)),
  winslash = "/",
  mustWork = TRUE
)

input_file <- file.path(root_dir, "inputs", "Galaxy452-gene_level_values.tabular")
analysis_dir <- file.path(root_dir, "edgeR_batch_corrected_analysis")
cache_dir <- file.path(analysis_dir, "cache")
table_dir <- file.path(analysis_dir, "tables")
figure_dir <- file.path(analysis_dir, "figures")
object_dir <- file.path(analysis_dir, "objects")
report_dir <- file.path(analysis_dir, "report")
log_dir <- file.path(analysis_dir, "logs")

for (d in c(analysis_dir, cache_dir, table_dir, figure_dir, object_dir, report_dir, log_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

force_run <- any(commandArgs(trailingOnly = TRUE) %in% c("--force", "-f")) ||
  identical(Sys.getenv("RNASEQ_FORCE"), "1")

stage_is_fresh <- function(target, deps) {
  if (force_run || !file.exists(target)) return(FALSE)
  deps <- deps[file.exists(deps)]
  if (length(deps) == 0) return(FALSE)
  file.info(target)$mtime >= max(file.info(deps)$mtime, na.rm = TRUE)
}

write_stage_stamp <- function(stage_name) {
  writeLines(
    c(
      paste("stage:", stage_name),
      paste("timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
      paste("root_dir:", root_dir)
    ),
    file.path(log_dir, paste0(stage_name, "_complete.txt"))
  )
}

read_stage <- function(path) {
  if (!file.exists(path)) stop("Required cached object is missing: ", path)
  readRDS(path)
}

save_gg <- function(plot, stem, width = 7, height = 5, dir = figure_dir) {
  ggsave(file.path(dir, paste0(stem, ".png")), plot, width = width, height = height,
         units = "in", dpi = 300, bg = "white")
  ggsave(file.path(dir, paste0(stem, ".pdf")), plot, width = width, height = height,
         units = "in", device = cairo_pdf, bg = "white")
}

base_theme <- theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold")
  )

condition_colors <- c(C = "#4D4D4D", `2DG` = "#D9822B", LDHA = "#2F6F9F")
condition_labels <- c(C = "Control", `2DG` = "2DG + OX", LDHA = "LDHA + OX")
deg_colors <- c(Down = "#2F6F9F", NS = "#BDBDBD", Up = "#D9822B", `FDR only` = "#7A7A7A")

contrast_specs <- list(
  `2DG_vs_C` = list(coef = "condition2DG", label = "2DG + OX vs Control", condition = "2DG"),
  `LDHA_vs_C` = list(coef = "conditionLDHA", label = "LDHA + OX vs Control", condition = "LDHA")
)

classify_deg <- function(res) {
  sig <- rep("NS", nrow(res))
  sig[res$FDR < 0.05 & abs(res$logFC) < 1] <- "FDR only"
  sig[res$FDR < 0.05 & res$logFC >= 1] <- "Up"
  sig[res$FDR < 0.05 & res$logFC <= -1] <- "Down"
  factor(sig, levels = c("Down", "NS", "FDR only", "Up"))
}

parse_lane_name <- function(x) {
  m <- regexec("^(PCS[0-9]+)_(C|2DG|LDHA)_(.+)$", x)
  hit <- regmatches(x, m)[[1]]
  if (length(hit) != 4) stop("Could not parse sample name: ", x)
  data.frame(
    lane_name = x,
    pcs = hit[2],
    condition = hit[3],
    technical_label = hit[4],
    stringsAsFactors = FALSE
  )
}

annotate_genes <- function(gene_ids) {
  ensembl_base <- sub("\\..*$", "", gene_ids)
  annotation <- data.frame(
    gene_id = gene_ids,
    ensembl_id = ensembl_base,
    stringsAsFactors = FALSE
  )
  ann_map <- suppressMessages(
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys = unique(ensembl_base),
      keytype = "ENSEMBL",
      columns = c("SYMBOL", "GENENAME", "ENTREZID")
    )
  )
  ann_map <- ann_map[order(is.na(ann_map$SYMBOL), is.na(ann_map$GENENAME)), ]
  ann_map <- ann_map[!duplicated(ann_map$ENSEMBL), ]
  annotation$SYMBOL <- ann_map$SYMBOL[match(annotation$ensembl_id, ann_map$ENSEMBL)]
  annotation$GENENAME <- ann_map$GENENAME[match(annotation$ensembl_id, ann_map$ENSEMBL)]
  annotation$ENTREZID <- ann_map$ENTREZID[match(annotation$ensembl_id, ann_map$ENSEMBL)]
  annotation$gene_label <- ifelse(!is.na(annotation$SYMBOL) & annotation$SYMBOL != "",
                                  annotation$SYMBOL, annotation$ensembl_id)
  annotation
}

message("RNA-seq pipeline root: ", root_dir)
