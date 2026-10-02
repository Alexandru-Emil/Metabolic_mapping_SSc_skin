# Scientific workflow derived from Revision/RNAseq PCS met inhibitors/results/scripts/03_functional_analysis.R
# Run from the repository root with METABOLIC_INPUT_DIR set to the external project data.
source("data_preprocessing/common/config.R")
script_dir <- {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)) else getwd()
}
source(file.path(script_dir, "00_config.R"))

stage_name <- "03_functional_analysis"
stage_rds <- file.path(cache_dir, paste0(stage_name, ".rds"))
stage1_rds <- file.path(cache_dir, "01_qc_batch_correction.rds")
stage2_rds <- file.path(cache_dir, "02_deg_edger.rds")
deps <- c(stage1_rds, stage2_rds, file.path(script_dir, "00_config.R"), file.path(script_dir, "03_functional_analysis.R"))

if (stage_is_fresh(stage_rds, deps)) {
  message("Using cached ", stage_name, ": ", stage_rds)
  quit(save = "no", status = 0)
}

message("Stage 3: GO Biological Process GSEA and ORA from edgeR results")
s1 <- read_stage(stage1_rds)
s2 <- read_stage(stage2_rds)
annotation <- s1$filtered_annotation
edgeR_results <- s2$edgeR_results
universe_entrez <- unique(na.omit(annotation$ENTREZID))

make_ranked_entrez <- function(res) {
  x <- res[!is.na(res$ENTREZID) & is.finite(res$rank_metric_signed_sqrt_F), c("ENTREZID", "rank_metric_signed_sqrt_F")]
  x <- x[order(abs(x$rank_metric_signed_sqrt_F), decreasing = TRUE), ]
  x <- x[!duplicated(x$ENTREZID), ]
  ranks <- x$rank_metric_signed_sqrt_F
  names(ranks) <- x$ENTREZID
  sort(ranks, decreasing = TRUE)
}

gsea_results <- list()
ora_results <- list()
for (nm in names(edgeR_results)) {
  res <- edgeR_results[[nm]]
  ranks <- make_ranked_entrez(res)
  gsea <- tryCatch({
    suppressMessages(
      gseGO(
        geneList = ranks,
        OrgDb = org.Hs.eg.db,
        ont = "BP",
        keyType = "ENTREZID",
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        verbose = FALSE,
        eps = 0
      )
    )
  }, error = function(e) {
    writeLines(conditionMessage(e), file.path(log_dir, paste0("edgeR_gseGO_error_", nm, ".txt")))
    NULL
  })
  gsea_df <- if (!is.null(gsea)) as.data.frame(gsea) else data.frame()
  gsea_results[[nm]] <- gsea_df
  write.csv(gsea_df, file.path(table_dir, paste0("edgeR_go_bp_gsea_", nm, ".csv")), row.names = FALSE)

  sig_up <- unique(na.omit(res$ENTREZID[res$FDR < 0.05 & res$logFC >= 1]))
  sig_down <- unique(na.omit(res$ENTREZID[res$FDR < 0.05 & res$logFC <= -1]))
  for (direction in c("up", "down")) {
    genes <- if (direction == "up") sig_up else sig_down
    ora_nm <- paste(nm, direction, sep = "_")
    if (length(genes) >= 10) {
      ora <- tryCatch({
        suppressMessages(
          enrichGO(
            gene = genes,
            universe = universe_entrez,
            OrgDb = org.Hs.eg.db,
            keyType = "ENTREZID",
            ont = "BP",
            pAdjustMethod = "BH",
            pvalueCutoff = 1,
            qvalueCutoff = 1,
            minGSSize = 10,
            maxGSSize = 500,
            readable = TRUE
          )
        )
      }, error = function(e) {
        writeLines(conditionMessage(e), file.path(log_dir, paste0("edgeR_enrichGO_error_", ora_nm, ".txt")))
        NULL
      })
      ora_df <- if (!is.null(ora)) as.data.frame(ora) else data.frame()
    } else {
      ora_df <- data.frame()
    }
    ora_results[[ora_nm]] <- ora_df
    write.csv(ora_df, file.path(table_dir, paste0("edgeR_go_bp_ora_", ora_nm, ".csv")), row.names = FALSE)
  }
}

plot_gsea <- function(df, stem, title) {
  if (nrow(df) == 0) {
    p <- ggplot() +
      annotate("text", x = 0, y = 0, label = "No GO BP GSEA results available", size = 4) +
      xlim(-1, 1) + ylim(-1, 1) +
      labs(title = title, x = NULL, y = NULL) +
      theme_void(base_size = 11)
    save_gg(p, stem, width = 7, height = 4.5)
    return(invisible(NULL))
  }
  df <- df[is.finite(df$NES) & !is.na(df$p.adjust), ]
  pos <- head(df[df$NES > 0, ][order(df$p.adjust[df$NES > 0], -abs(df$NES[df$NES > 0])), ], 8)
  neg <- head(df[df$NES < 0, ][order(df$p.adjust[df$NES < 0], -abs(df$NES[df$NES < 0])), ], 8)
  top <- rbind(pos, neg)
  if (nrow(top) == 0) top <- head(df[order(df$p.adjust), ], 16)
  top$Description <- factor(top$Description, levels = top$Description[order(top$NES)])
  p <- ggplot(top, aes(x = NES, y = Description, size = setSize, color = p.adjust)) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.35) +
    geom_point(alpha = 0.9) +
    scale_color_gradient(low = "#D9822B", high = "#2F6F9F", trans = "reverse") +
    labs(
      title = title,
      subtitle = "Top positive and negative GO Biological Process GSEA terms",
      x = "Normalized enrichment score",
      y = NULL,
      color = "Adjusted P",
      size = "Gene set size"
    ) +
    base_theme
  save_gg(p, stem, width = 8.5, height = 6.3)
}

for (nm in names(gsea_results)) {
  plot_gsea(gsea_results[[nm]], paste0("pathway_edgeR_go_bp_gsea_", nm),
            paste0("GO BP GSEA: ", contrast_specs[[nm]]$label))
}

top_gsea_summary <- do.call(rbind, lapply(names(gsea_results), function(nm) {
  df <- gsea_results[[nm]]
  if (nrow(df) == 0) return(data.frame())
  df <- df[order(df$p.adjust, -abs(df$NES)), ]
  df <- head(df[, c("ID", "Description", "setSize", "NES", "pvalue", "p.adjust", "qvalue")], 10)
  df$contrast <- nm
  df$contrast_label <- contrast_specs[[nm]]$label
  df
}))
write.csv(top_gsea_summary, file.path(table_dir, "edgeR_top_go_gsea_summary_for_report.csv"), row.names = FALSE)

stage_data <- list(
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  pathway_model = "clusterProfiler GO Biological Process GSEA on edgeR signed sqrt(QL F-statistic); ORA on high-confidence DEG sets",
  gsea_results = gsea_results,
  ora_results = ora_results,
  top_gsea_summary = top_gsea_summary
)
saveRDS(stage_data, stage_rds)
write_stage_stamp(stage_name)
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_03_functional_analysis.txt"))
message("Stage 3 complete: ", stage_rds)
