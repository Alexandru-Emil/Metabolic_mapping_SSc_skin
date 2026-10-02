# Scientific workflow derived from Revision/RNAseq PCS met inhibitors/results/scripts/06_combatseq_deg_go_reports.R
# Run from the repository root with METABOLIC_INPUT_DIR set to the external project data.
source("data_preprocessing/common/config.R")
script_dir <- {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)) else getwd()
}
source(file.path(script_dir, "00_config.R"))

stage_name <- "06_combatseq_deg_go_reports"
stage_rds <- file.path(cache_dir, paste0(stage_name, ".rds"))
stage1_rds <- file.path(cache_dir, "01_qc_batch_correction.rds")
stage5_rds <- file.path(cache_dir, "05_method_threshold_sensitivity.rds")
deps <- c(
  stage1_rds,
  stage5_rds,
  file.path(script_dir, "00_config.R"),
  file.path(script_dir, "06_combatseq_deg_go_reports.R")
)

if (stage_is_fresh(stage_rds, deps)) {
  message("Using cached ", stage_name, ": ", stage_rds)
  quit(save = "no", status = 0)
}

message("Stage 6: method-specific ComBat-seq DEG and GO reports")
s1 <- read_stage(stage1_rds)
s5 <- read_stage(stage5_rds)

fc_threshold <- 1.5
abs_log2fc_threshold <- log2(fc_threshold)
fdr_threshold <- 0.05
go_onts <- c("BP", "MF", "CC")
go_directions <- c("All", "Up", "Down")

method_specs <- list(
  edgeR_exact_ComBatSeq = list(
    method = "edgeR_exact_ComBatSeq",
    label = "edgeR exact test on ComBat-seq counts",
    short = "edgeR exact ComBat-seq",
    prefix = "edgeR_exact_ComBatSeq_FDR05_FC15",
    report_title = "edgeR exact ComBat-seq DEG and GO analysis",
    method_note = "edgeR exactTest was run on condition-preserving ComBat-seq PCS-corrected counts, comparing each treatment with control."
  ),
  DESeq2_ComBatSeq = list(
    method = "DESeq2_ComBatSeq",
    label = "DESeq2 on ComBat-seq counts",
    short = "DESeq2 ComBat-seq",
    prefix = "DESeq2_ComBatSeq_FDR05_FC15",
    report_title = "DESeq2 ComBat-seq DEG and GO analysis",
    method_note = "DESeq2 Wald tests were run on condition-preserving ComBat-seq PCS-corrected counts with design ~ condition."
  )
)

sample_meta <- s1$sample_meta
sample_meta <- sample_meta[sample_meta$sample_id, , drop = FALSE]
universe_entrez <- unique(na.omit(s1$filtered_annotation$ENTREZID))

wrap_label <- function(x, width = 58) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))
}

ratio_to_numeric <- function(x) {
  vapply(strsplit(as.character(x), "/", fixed = TRUE), function(parts) {
    if (length(parts) != 2) return(NA_real_)
    as.numeric(parts[1]) / as.numeric(parts[2])
  }, numeric(1))
}

safe_sheet <- local({
  used <- character()
  function(name) {
    x <- gsub("[\\[\\]\\*\\?/\\\\:]", "_", name)
    x <- substr(x, 1, 31)
    base <- x
    i <- 1
    while (tolower(x) %in% tolower(used)) {
      suffix <- paste0("_", i)
      x <- paste0(substr(base, 1, 31 - nchar(suffix)), suffix)
      i <- i + 1
    }
    used <<- c(used, x)
    x
  }
})

write_table_sheet <- function(wb, sheet, df) {
  sheet <- safe_sheet(sheet)
  addWorksheet(wb, sheet)
  if (is.null(df) || nrow(df) == 0) {
    writeData(wb, sheet, data.frame(note = "No rows returned for this table."))
  } else {
    writeDataTable(wb, sheet, df)
  }
}

classify_requested_deg <- function(res) {
  out <- rep("NS", nrow(res))
  out[res$FDR < fdr_threshold & res$logFC >= abs_log2fc_threshold] <- "Up"
  out[res$FDR < fdr_threshold & res$logFC <= -abs_log2fc_threshold] <- "Down"
  factor(out, levels = c("Down", "NS", "Up"))
}

make_deg_summary <- function(method_df) {
  do.call(rbind, lapply(names(contrast_specs), function(nm) {
    res <- method_df[method_df$contrast == nm, ]
    deg <- res[res$is_deg, ]
    top <- res[order(res$FDR, res$PValue, -abs(res$logFC), na.last = TRUE), ][1, ]
    data.frame(
      contrast = nm,
      contrast_label = contrast_specs[[nm]]$label,
      genes_tested = nrow(res),
      up = sum(deg$direction == "Up", na.rm = TRUE),
      down = sum(deg$direction == "Down", na.rm = TRUE),
      total_deg = nrow(deg),
      min_pvalue = min(res$PValue, na.rm = TRUE),
      min_fdr = min(res$FDR, na.rm = TRUE),
      top_gene = top$gene_label,
      top_gene_logFC = top$logFC,
      top_gene_PValue = top$PValue,
      top_gene_FDR = top$FDR,
      stringsAsFactors = FALSE
    )
  }))
}

run_go_ora <- function(genes_entrez, universe_entrez, ont, id) {
  genes_entrez <- unique(na.omit(genes_entrez))
  if (length(genes_entrez) < 10) return(data.frame())
  out <- tryCatch({
    suppressMessages(
      enrichGO(
        gene = genes_entrez,
        universe = universe_entrez,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = ont,
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        minGSSize = 10,
        maxGSSize = 500,
        readable = TRUE
      )
    )
  }, error = function(e) {
    writeLines(conditionMessage(e), file.path(log_dir, paste0("GO_ORA_error_", id, "_", ont, ".txt")))
    NULL
  })
  if (is.null(out)) data.frame() else as.data.frame(out)
}

plot_deg_bar <- function(summary_df, stem, title) {
  plot_df <- rbind(
    data.frame(contrast_label = summary_df$contrast_label, direction = "Up", n = summary_df$up),
    data.frame(contrast_label = summary_df$contrast_label, direction = "Down", n = summary_df$down)
  )
  plot_df$direction <- factor(plot_df$direction, levels = c("Down", "Up"))
  p <- ggplot(plot_df, aes(x = contrast_label, y = n, fill = direction)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.62, color = "grey25", linewidth = 0.25) +
    geom_text(aes(label = n), position = position_dodge(width = 0.72), vjust = -0.35, size = 3.2) +
    scale_fill_manual(values = c(Down = "#2F6F9F", Up = "#D9822B")) +
    labs(
      title = title,
      subtitle = "DEGs use FDR < 0.05 and |log2FC| >= log2(1.5)",
      x = NULL,
      y = "Number of DEGs",
      fill = "Direction"
    ) +
    coord_cartesian(ylim = c(0, max(plot_df$n, 1) * 1.16)) +
    base_theme
  save_gg(p, stem, width = 7.4, height = 4.7)
}

plot_volcano <- function(res, stem, title) {
  res$neg_log10_fdr <- -log10(pmax(res$FDR, .Machine$double.xmin))
  labels <- res[res$is_deg, ]
  labels <- labels[order(labels$FDR, -abs(labels$logFC)), ]
  labels <- head(labels, 18)
  p <- ggplot(res, aes(x = logFC, y = neg_log10_fdr, color = direction)) +
    geom_point(size = 1.05, alpha = 0.72, stroke = 0) +
    geom_vline(xintercept = c(-abs_log2fc_threshold, abs_log2fc_threshold),
               linetype = "dashed", color = "grey45", linewidth = 0.35) +
    geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed", color = "grey45", linewidth = 0.35) +
    ggrepel::geom_text_repel(
      data = labels,
      aes(label = gene_label),
      size = 2.8,
      max.overlaps = Inf,
      box.padding = 0.34,
      min.segment.length = 0,
      show.legend = FALSE
    ) +
    scale_color_manual(values = c(Down = "#2F6F9F", NS = "#BDBDBD", Up = "#D9822B"), drop = FALSE) +
    labs(
      title = title,
      subtitle = "Dashed lines mark FDR 0.05 and FC 1.5",
      x = "log2 fold-change",
      y = "-log10 FDR",
      color = "Class"
    ) +
    base_theme
  save_gg(p, stem, width = 7.3, height = 5.7)
}

plot_go_dot <- function(go_df, stem, title) {
  bp <- go_df[go_df$ont == "BP" & go_df$direction %in% go_directions, ]
  bp <- bp[is.finite(bp$p.adjust) & !is.na(bp$p.adjust), ]
  if (nrow(bp) == 0) {
    p <- ggplot() +
      annotate("text", x = 0, y = 0, label = "No GO BP enrichment terms returned", size = 4) +
      xlim(-1, 1) + ylim(-1, 1) +
      labs(title = title, x = NULL, y = NULL) +
      theme_void(base_size = 11)
    save_gg(p, stem, width = 8, height = 4.4)
    return(invisible(NULL))
  }
  top <- do.call(rbind, lapply(split(bp, bp$direction), function(df) {
    df <- df[order(df$p.adjust, -df$Count), ]
    head(df, 8)
  }))
  top$gene_ratio_numeric <- ratio_to_numeric(top$GeneRatio)
  top$Description_wrapped <- wrap_label(top$Description)
  top$Description_wrapped <- factor(top$Description_wrapped, levels = rev(unique(top$Description_wrapped)))
  top$direction <- factor(top$direction, levels = go_directions)
  p <- ggplot(top, aes(x = gene_ratio_numeric, y = Description_wrapped, size = Count, color = p.adjust)) +
    geom_point(alpha = 0.9) +
    facet_grid(direction ~ ., scales = "free_y", space = "free_y") +
    scale_color_gradient(low = "#D9822B", high = "#2F6F9F", trans = "reverse") +
    labs(
      title = title,
      subtitle = "Top GO Biological Process terms from DEG over-representation analysis",
      x = "Gene ratio",
      y = NULL,
      color = "Adjusted P",
      size = "Genes"
    ) +
    base_theme +
    theme(strip.text.y = element_text(angle = 0))
  save_gg(p, stem, width = 8.8, height = max(5.4, 1.4 + 0.34 * nrow(top)))
}

plot_deg_heatmap <- function(method_df, stem, title) {
  degs <- method_df[method_df$is_deg, ]
  top_ids <- unique(unlist(lapply(split(degs, degs$contrast), function(df) {
    df <- df[order(df$FDR, -abs(df$logFC)), ]
    head(df$gene_id, 35)
  })))
  if (length(top_ids) < 2) {
    p <- ggplot() +
      annotate("text", x = 0, y = 0, label = "Fewer than two DEGs available for heatmap", size = 4) +
      xlim(-1, 1) + ylim(-1, 1) +
      labs(title = title, x = NULL, y = NULL) +
      theme_void(base_size = 11)
    save_gg(p, stem, width = 7, height = 4.5)
    return(invisible(NULL))
  }
  top_ids <- intersect(top_ids, rownames(s1$batch_corrected_logcpm))
  heat_mat <- s1$batch_corrected_logcpm[top_ids, sample_meta$sample_id, drop = FALSE]
  labels <- method_df$gene_label[match(rownames(heat_mat), method_df$gene_id)]
  labels <- make.unique(ifelse(is.na(labels) | labels == "", rownames(heat_mat), labels))
  rownames(heat_mat) <- labels
  ann_col <- data.frame(PCS = sample_meta$pcs, Condition = sample_meta$condition)
  rownames(ann_col) <- sample_meta$sample_id
  ann_colors <- list(
    Condition = condition_colors,
    PCS = setNames(brewer.pal(max(3, length(levels(sample_meta$pcs))), "Set2")[seq_along(levels(sample_meta$pcs))],
                   levels(sample_meta$pcs))
  )
  pheatmap(
    heat_mat,
    scale = "row",
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    color = colorRampPalette(c("#2166ac", "#f7f7f7", "#b35806"))(100),
    border_color = NA,
    fontsize_row = ifelse(nrow(heat_mat) > 70, 5.5, 7),
    fontsize_col = 8,
    main = title,
    filename = file.path(figure_dir, paste0(stem, ".png")),
    width = 8.6,
    height = max(6.6, min(13, 2 + nrow(heat_mat) * 0.13))
  )
  pheatmap(
    heat_mat,
    scale = "row",
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    color = colorRampPalette(c("#2166ac", "#f7f7f7", "#b35806"))(100),
    border_color = NA,
    fontsize_row = ifelse(nrow(heat_mat) > 70, 5.5, 7),
    fontsize_col = 8,
    main = title,
    filename = file.path(figure_dir, paste0(stem, ".pdf")),
    width = 8.6,
    height = max(6.6, min(13, 2 + nrow(heat_mat) * 0.13))
  )
}

build_method_report <- function(spec, method_df, deg_summary, top_degs, go_results, go_summary, report_data_path) {
  rmd_path <- file.path(report_dir, paste0(spec$prefix, "_report.Rmd"))
  report_lines <- c(
    "---",
    paste0("title: \"", spec$report_title, "\""),
    "output:",
    "  html_document:",
    "    toc: true",
    "    toc_depth: 3",
    "    number_sections: false",
    "    theme: cosmo",
    "    df_print: paged",
    "params:",
    "  data_path: \"\"",
    "  out_dir: \"\"",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
    "report <- readRDS(params$data_path)",
    "out_dir <- params$out_dir",
    "fig <- function(name) knitr::include_graphics(file.path(out_dir, 'figures', name))",
    "fmt_int <- function(x) format(x, big.mark = ',', scientific = FALSE, trim = TRUE)",
    "deg_summary <- report$deg_summary",
    "top_degs <- report$top_degs",
    "go_summary <- report$go_summary",
    "```",
    "",
    "## Technical summary",
    "",
    paste0("This report uses **", spec$label, "** as requested. DEGs are defined as genes with FDR < 0.05 and absolute fold-change >= 1.5, equivalent to |log2FC| >= ", round(abs_log2fc_threshold, 4), ". GO results use over-representation analysis on the resulting DEG sets, with the filtered RNA-seq gene universe as background."),
    "",
    paste0(spec$method_note, " The analysis uses the cached ComBat-seq matrix created with PCS as batch and condition preserved."),
    "",
    "```{r deg-summary}",
    "show <- deg_summary[, c('contrast_label', 'genes_tested', 'up', 'down', 'total_deg', 'min_fdr', 'top_gene', 'top_gene_FDR')]",
    "colnames(show) <- c('Contrast', 'Genes tested', 'Up', 'Down', 'Total DEG', 'Minimum FDR', 'Top gene', 'Top gene FDR')",
    "knitr::kable(show, digits = 4)",
    "```",
    "",
    "## DEG burden and direction are threshold-dependent",
    "",
    "The bar chart counts only genes passing the requested threshold. Direction is based on the sign of the treatment-versus-control log2 fold-change.",
    "",
    "```{r deg-bar, out.width='100%'}",
    "fig(report$figures$deg_bar)",
    "```",
    "",
    "## Volcano plots show the selected DEG boundary",
    "",
    "These volcano plots place the requested FDR and fold-change cutoffs directly on the all-gene result. Genes labeled in the plot are the strongest DEG calls by FDR and effect size.",
    "",
    "### 2DG + OX versus control",
    "",
    "```{r volcano-2dg, out.width='100%'}",
    "fig(report$figures$volcano[['2DG_vs_C']])",
    "```",
    "",
    "### LDHA + OX versus control",
    "",
    "```{r volcano-ldha, out.width='100%'}",
    "fig(report$figures$volcano[['LDHA_vs_C']])",
    "```",
    "",
    "## PCS-corrected expression separates top DEG patterns",
    "",
    "The heatmap uses ComBat-seq PCS-corrected logCPM for the top DEG set across both contrasts. It is an expression-display figure; exact statistics should be read from the all-gene and DEG tables.",
    "",
    "```{r heatmap, out.width='100%'}",
    "fig(report$figures$heatmap)",
    "```",
    "",
    "```{r top-degs}",
    "td <- top_degs[, c('contrast_label', 'direction', 'gene_label', 'GENENAME', 'logFC', 'PValue', 'FDR')]",
    "colnames(td) <- c('Contrast', 'Direction', 'Gene', 'Gene name', 'log2FC', 'P value', 'FDR')",
    "knitr::kable(td, digits = 4)",
    "```",
    "",
    "## GO enrichment uses the selected DEG sets",
    "",
    "GO over-representation analysis was run separately for all, up-regulated, and down-regulated DEG sets. The plotted terms show GO Biological Process results; full BP, MF, and CC tables are included in the workbook and CSV outputs.",
    "",
    "### 2DG + OX versus control",
    "",
    "```{r go-2dg, out.width='100%'}",
    "fig(report$figures$go_bp[['2DG_vs_C']])",
    "```",
    "",
    "### LDHA + OX versus control",
    "",
    "```{r go-ldha, out.width='100%'}",
    "fig(report$figures$go_bp[['LDHA_vs_C']])",
    "```",
    "",
    "```{r go-summary}",
    "if (nrow(go_summary) > 0) {",
    "  gs <- go_summary[, c('contrast_label', 'direction', 'ont', 'Description', 'GeneRatio', 'Count', 'p.adjust', 'qvalue')]",
    "  colnames(gs) <- c('Contrast', 'DEG set', 'Ontology', 'GO term', 'Gene ratio', 'Genes', 'FDR', 'q value')",
    "  knitr::kable(gs, digits = 4)",
    "} else {",
    "  cat('No GO terms were returned for the selected DEG sets.')",
    "}",
    "```",
    "",
    "## Scope, data, and metric definitions",
    "",
    "- Comparisons: 2DG + OX vs control and LDHA + OX vs control.",
    "- Batch correction: ComBat-seq was run with PCS as the batch variable and condition preserved.",
    "- DEG threshold: FDR < 0.05 and FC >= 1.5, implemented as |log2FC| >= log2(1.5).",
    "- GO universe: all retained, filtered RNA-seq genes that could be mapped to Entrez identifiers.",
    "- Direction: positive log2FC indicates higher expression in the treatment condition relative to control.",
    "",
    "## Methodology and model contract",
    "",
    paste0("- ", spec$method_note),
    "- DEG tables were filtered after model fitting; all-gene results are exported so alternative thresholds remain auditable.",
    "- GO over-representation analysis used `clusterProfiler::enrichGO` with Benjamini-Hochberg adjustment and the filtered expressed gene universe.",
    "- Reports and tables are generated from cached RDS objects so reruns skip unchanged upstream QC and DEG stages.",
    "",
    "## Limitations and interpretation",
    "",
    "- These are method-specific reports requested on ComBat-seq corrected counts. They should not be mixed without qualification with the primary edgeR QL model that includes PCS directly in the negative-binomial design.",
    "- Correcting counts before DEG testing can alter variance structure and may make corrected-count screens more permissive than a paired count model.",
    "- After ComBat-seq correction, the pair/block structure is not directly modeled in the exact-test or condition-only DESeq2 comparison.",
    "- GO over-representation results depend on DEG threshold, annotation coverage, and redundancy among GO terms.",
    "- Upstream FASTQ/BAM QC, RNA quality, alignment rates, strandedness, and library-prep covariates were not supplied.",
    "",
    "## Recommended next steps",
    "",
    "- Use the method-specific DEG tables in this report when the manuscript or supplement explicitly refers to ComBat-seq corrected-count DEG analysis.",
    "- Prioritize genes and GO terms that are robust across the edgeR exact and DESeq2 ComBat-seq reports.",
    "- Keep the previous edgeR QL PCS-adjusted analysis as the conservative primary inferential reference, and present these results as requested sensitivity or corrected-count reports."
  )
  writeLines(report_lines, rmd_path)
  Sys.setenv(RSTUDIO_PANDOC = project_path("external/tools"))
  rmarkdown::render(
    input = rmd_path,
    output_file = paste0(spec$prefix, "_report.html"),
    output_dir = report_dir,
    params = list(data_path = report_data_path, out_dir = analysis_dir),
    quiet = TRUE,
    envir = new.env(parent = globalenv())
  )
  file.path(report_dir, paste0(spec$prefix, "_report.html"))
}

method_outputs <- list()

for (method_key in names(method_specs)) {
  spec <- method_specs[[method_key]]
  message("Building ", spec$label)
  method_df <- s5$all_method_results[as.character(s5$all_method_results$method) == spec$method, ]
  if (nrow(method_df) == 0) stop("No cached method results found for ", spec$method)
  method_df$is_deg <- !is.na(method_df$FDR) & method_df$FDR < fdr_threshold & abs(method_df$logFC) >= abs_log2fc_threshold
  method_df$direction <- classify_requested_deg(method_df)
  method_df$direction <- as.character(method_df$direction)

  deg_summary <- make_deg_summary(method_df)
  top_degs <- do.call(rbind, lapply(names(contrast_specs), function(nm) {
    deg <- method_df[method_df$contrast == nm & method_df$is_deg, ]
    if (nrow(deg) == 0) return(data.frame())
    deg <- deg[order(deg$FDR, deg$PValue, -abs(deg$logFC), na.last = TRUE), ]
    head(deg, 30)
  }))

  write.csv(method_df, file.path(table_dir, paste0(spec$prefix, "_all_genes.csv")), row.names = FALSE)
  write.csv(method_df[method_df$is_deg, ], file.path(table_dir, paste0(spec$prefix, "_DEGs_all_contrasts.csv")), row.names = FALSE)
  write.csv(deg_summary, file.path(table_dir, paste0(spec$prefix, "_DEG_summary.csv")), row.names = FALSE)
  write.csv(top_degs, file.path(table_dir, paste0(spec$prefix, "_top_DEGs_for_report.csv")), row.names = FALSE)

  plot_deg_bar(deg_summary, paste0(spec$prefix, "_DEG_summary_bar"), paste0(spec$short, ": DEG burden"))
  plot_deg_heatmap(method_df, paste0(spec$prefix, "_top_DEG_heatmap"), paste0(spec$short, ": top DEGs on PCS-corrected logCPM"))

  volcano_files <- list()
  go_plot_files <- list()
  go_results_list <- list()
  for (nm in names(contrast_specs)) {
    res <- method_df[method_df$contrast == nm, ]
    volcano_stem <- paste0(spec$prefix, "_volcano_", nm)
    plot_volcano(res, volcano_stem, paste0(spec$short, ": ", contrast_specs[[nm]]$label))
    volcano_files[[nm]] <- paste0(volcano_stem, ".png")

    deg <- res[res$is_deg, ]
    for (direction in go_directions) {
      if (direction == "All") {
        genes <- deg$ENTREZID
      } else {
        genes <- deg$ENTREZID[deg$direction == direction]
      }
      for (ont in go_onts) {
        go_id <- paste(spec$prefix, nm, direction, ont, sep = "_")
        go_df <- run_go_ora(genes, universe_entrez, ont, go_id)
        if (nrow(go_df) > 0) {
          go_df$method <- spec$method
          go_df$method_label <- spec$label
          go_df$contrast <- nm
          go_df$contrast_label <- contrast_specs[[nm]]$label
          go_df$direction <- direction
          go_df$ont <- ont
        }
        go_results_list[[go_id]] <- go_df
        write.csv(go_df, file.path(table_dir, paste0(go_id, "_GO_ORA.csv")), row.names = FALSE)
      }
    }

    go_combined_for_plot <- do.call(rbind, go_results_list[grepl(paste0("^", spec$prefix, "_", nm, "_"), names(go_results_list))])
    go_stem <- paste0(spec$prefix, "_GO_BP_ORA_", nm)
    plot_go_dot(go_combined_for_plot, go_stem, paste0(spec$short, ": GO BP ORA for ", contrast_specs[[nm]]$label))
    go_plot_files[[nm]] <- paste0(go_stem, ".png")
  }

  go_results <- do.call(rbind, go_results_list)
  if (is.null(go_results)) go_results <- data.frame()
  write.csv(go_results, file.path(table_dir, paste0(spec$prefix, "_GO_ORA_all_results.csv")), row.names = FALSE)
  go_summary <- if (nrow(go_results) > 0) {
    do.call(rbind, lapply(split(go_results, list(go_results$contrast, go_results$direction, go_results$ont), drop = TRUE), function(df) {
      df <- df[order(df$p.adjust, -df$Count), ]
      head(df, 5)
    }))
  } else {
    data.frame()
  }
  write.csv(go_summary, file.path(table_dir, paste0(spec$prefix, "_GO_ORA_top_terms_for_report.csv")), row.names = FALSE)

  wb <- createWorkbook()
  write_table_sheet(wb, "DEG_summary", deg_summary)
  write_table_sheet(wb, "top_DEGs", top_degs)
  for (nm in names(contrast_specs)) {
    res <- method_df[method_df$contrast == nm, ]
    write_table_sheet(wb, paste0("all_", nm), res)
    write_table_sheet(wb, paste0("DEG_", nm), res[res$is_deg, ])
    write_table_sheet(wb, paste0("DEG_up_", nm), res[res$is_deg & res$direction == "Up", ])
    write_table_sheet(wb, paste0("DEG_down_", nm), res[res$is_deg & res$direction == "Down", ])
  }
  write_table_sheet(wb, "GO_top_terms", go_summary)
  for (nm in names(go_results_list)) {
    write_table_sheet(wb, substr(nm, nchar(spec$prefix) + 2, nchar(nm)), go_results_list[[nm]])
  }
  workbook_path <- file.path(table_dir, paste0(spec$prefix, "_DEG_GO_results.xlsx"))
  saveWorkbook(wb, workbook_path, overwrite = TRUE)

  report_data <- list(
    method_spec = spec,
    fdr_threshold = fdr_threshold,
    fc_threshold = fc_threshold,
    abs_log2fc_threshold = abs_log2fc_threshold,
    deg_summary = deg_summary,
    top_degs = top_degs,
    go_summary = go_summary,
    figures = list(
      deg_bar = paste0(spec$prefix, "_DEG_summary_bar.png"),
      heatmap = paste0(spec$prefix, "_top_DEG_heatmap.png"),
      volcano = volcano_files,
      go_bp = go_plot_files
    ),
    workbook_path = workbook_path
  )
  report_data_path <- file.path(object_dir, paste0(spec$prefix, "_report_data.rds"))
  saveRDS(report_data, report_data_path)
  report_html <- build_method_report(spec, method_df, deg_summary, top_degs, go_results, go_summary, report_data_path)

  method_outputs[[method_key]] <- list(
    report_html = report_html,
    workbook_path = workbook_path,
    deg_summary = deg_summary,
    go_terms = nrow(go_results),
    report_data_path = report_data_path
  )
}

stage_data <- list(
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  threshold = list(FDR = fdr_threshold, fold_change = fc_threshold, abs_log2FC = abs_log2fc_threshold),
  method_outputs = method_outputs
)
saveRDS(stage_data, stage_rds)
write_stage_stamp(stage_name)
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_06_combatseq_deg_go_reports.txt"))
message("Stage 6 complete: ", stage_rds)
