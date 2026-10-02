# Scientific workflow derived from Revision/RNAseq PCS met inhibitors/results/scripts/07_symbol_interactive_reports.R
# Run from the repository root with METABOLIC_INPUT_DIR set to the external project data.
source("data_preprocessing/common/config.R")
script_dir <- {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)) else getwd()
}
source(file.path(script_dir, "00_config.R"))

stage_name <- "07_symbol_interactive_reports"
stage_rds <- file.path(cache_dir, paste0(stage_name, ".rds"))
stage1_rds <- file.path(cache_dir, "01_qc_batch_correction.rds")
stage5_rds <- file.path(cache_dir, "05_method_threshold_sensitivity.rds")
stage6_rds <- file.path(cache_dir, "06_combatseq_deg_go_reports.rds")
app_template <- file.path(script_dir, "symbol_report_app.R")

fdr_threshold <- 0.05
fc_threshold <- 1.5
abs_log2fc_threshold <- log2(fc_threshold)
version_tag <- "symbol_interactive_v2"

method_specs <- list(
  edgeR_exact_ComBatSeq = list(
    method = "edgeR_exact_ComBatSeq",
    title = "edgeR exact ComBat-seq symbol report",
    subtitle = "Differential expression and Gene Ontology results for 2DG + OX and LDHA + OX versus control.",
    prefix = "edgeR_exact_ComBatSeq_FDR05_FC15",
    method_label = "edgeR exact test on ComBat-seq corrected counts",
    method_note = "edgeR exactTest was run on condition-preserving ComBat-seq PCS-corrected counts as a two-group comparison for each treatment versus control."
  ),
  DESeq2_ComBatSeq = list(
    method = "DESeq2_ComBatSeq",
    title = "DESeq2 ComBat-seq symbol report",
    subtitle = "Differential expression and Gene Ontology results for 2DG + OX and LDHA + OX versus control.",
    prefix = "DESeq2_ComBatSeq_FDR05_FC15",
    method_label = "DESeq2 Wald test on ComBat-seq corrected counts",
    method_note = "DESeq2 Wald tests were run on condition-preserving ComBat-seq PCS-corrected counts with design ~ condition for each treatment versus control."
  )
)

app_dirs <- vapply(method_specs, function(spec) {
  file.path(report_dir, paste0(spec$prefix, "_", version_tag))
}, character(1))
expected_outputs <- unlist(lapply(app_dirs, function(path) {
  file.path(path, c("app.R", "report_data.rds", "start_report.R", "report_manifest.txt"))
}))

deps <- c(
  stage1_rds,
  stage5_rds,
  stage6_rds,
  app_template,
  file.path(script_dir, "00_config.R"),
  file.path(script_dir, "07_symbol_interactive_reports.R")
)

if (stage_is_fresh(stage_rds, deps) && all(file.exists(expected_outputs))) {
  message("Using cached ", stage_name, ": ", stage_rds)
  quit(save = "no", status = 0)
}

message("Stage 7: symbol-only interactive DEG and GO reports")
s1 <- read_stage(stage1_rds)
s5 <- read_stage(stage5_rds)

fill_missing_symbols <- function(data) {
  missing <- is.na(data$SYMBOL) | trimws(data$SYMBOL) == ""
  if (any(missing)) {
    keys_to_map <- unique(data$ensembl_id[missing])
    keys_to_map <- keys_to_map[!is.na(keys_to_map) & nzchar(keys_to_map)]
    mapped <- tryCatch(
      suppressMessages(
        AnnotationDbi::mapIds(
          org.Hs.eg.db,
          keys = keys_to_map,
          column = "SYMBOL",
          keytype = "ENSEMBL",
          multiVals = "first"
        )
      ),
      error = function(e) NULL
    )
    if (!is.null(mapped)) {
      data$SYMBOL[missing] <- unname(mapped[data$ensembl_id[missing]])
    }
  }
  data
}

clean_gene_name <- function(x) {
  x[is.na(x)] <- ""
  trimws(x)
}

clean_symbol_list <- function(x) {
  vapply(strsplit(as.character(x), "/", fixed = TRUE), function(tokens) {
    tokens <- trimws(tokens)
    tokens <- tokens[nzchar(tokens) & !grepl("^ENSG[0-9]+", tokens, ignore.case = TRUE)]
    paste(unique(tokens), collapse = ", ")
  }, character(1))
}

make_summary <- function(method_data) {
  do.call(rbind, lapply(names(contrast_specs), function(contrast) {
    all_rows <- method_data[method_data$contrast == contrast, , drop = FALSE]
    called <- all_rows[all_rows$is_DEG, , drop = FALSE]
    mapped <- called[called$has_symbol, , drop = FALSE]
    mapped <- mapped[order(mapped$FDR, -abs(mapped$logFC)), , drop = FALSE]
    data.frame(
      Contrast_key = contrast,
      Contrast = contrast_specs[[contrast]]$label,
      Genes_tested = nrow(all_rows),
      Called_DEGs = nrow(called),
      Symbol_mapped_DEGs = nrow(mapped),
      Unmapped_DEGs = nrow(called) - nrow(mapped),
      Mapped_Up = sum(mapped$Direction == "Up"),
      Mapped_Down = sum(mapped$Direction == "Down"),
      Top_symbol = if (nrow(mapped) > 0) mapped$SYMBOL[[1]] else "None",
      stringsAsFactors = FALSE
    )
  }))
}

make_symbol_gene_tables <- function(method_data) {
  mapped <- method_data[method_data$has_symbol, , drop = FALSE]
  mapped <- mapped[order(
    mapped$contrast,
    mapped$SYMBOL,
    is.na(mapped$FDR),
    mapped$FDR,
    is.na(mapped$PValue),
    mapped$PValue,
    -abs(mapped$logFC)
  ), , drop = FALSE]
  mapped <- mapped[!duplicated(paste(mapped$contrast, mapped$SYMBOL, sep = "||")), , drop = FALSE]
  mapped$row_key <- seq_len(nrow(mapped))

  all_genes <- data.frame(
    row_key = mapped$row_key,
    Symbol = mapped$SYMBOL,
    Gene_name = clean_gene_name(mapped$GENENAME),
    Contrast_key = mapped$contrast,
    Contrast = mapped$contrast_label,
    Direction = mapped$Direction,
    log2FC = mapped$logFC,
    Fold_change = 2^abs(mapped$logFC),
    P_value = mapped$PValue,
    FDR = mapped$FDR,
    Average_expression = mapped$average_expression,
    is_DEG = mapped$is_DEG,
    stringsAsFactors = FALSE
  )

  degs <- all_genes[all_genes$is_DEG, , drop = FALSE]
  degs <- degs[order(degs$Contrast_key, degs$FDR, -abs(degs$log2FC)), , drop = FALSE]

  expression_rows <- mapped[mapped$is_DEG, c("row_key", "gene_id", "SYMBOL"), drop = FALSE]
  expression_rows <- expression_rows[expression_rows$gene_id %in% rownames(s1$batch_corrected_logcpm), , drop = FALSE]
  expression <- if (nrow(expression_rows) > 0) {
    do.call(rbind, lapply(seq_len(nrow(expression_rows)), function(i) {
      gene_id <- expression_rows$gene_id[[i]]
      data.frame(
        row_key = expression_rows$row_key[[i]],
        Symbol = expression_rows$SYMBOL[[i]],
        Sample = colnames(s1$batch_corrected_logcpm),
        logCPM = as.numeric(s1$batch_corrected_logcpm[gene_id, ]),
        stringsAsFactors = FALSE
      )
    }))
  } else {
    data.frame(row_key = integer(), Symbol = character(), Sample = character(), logCPM = numeric())
  }

  list(all_genes = all_genes, degs = degs, expression = expression)
}

make_go_table <- function(path) {
  raw <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(raw) == 0) {
    return(data.frame(
      GO_ID = character(), Term = character(), Ontology = character(), DEG_set = character(),
      Contrast_key = character(), Contrast = character(), Gene_ratio = character(),
      Background_ratio = character(), Fold_enrichment = numeric(), Gene_count = integer(),
      P_value = numeric(), FDR = numeric(), q_value = numeric(), Genes = character()
    ))
  }
  data.frame(
    GO_ID = raw$ID,
    Term = raw$Description,
    Ontology = raw$ont,
    DEG_set = raw$direction,
    Contrast_key = raw$contrast,
    Contrast = raw$contrast_label,
    Gene_ratio = raw$GeneRatio,
    Background_ratio = raw$BgRatio,
    Fold_enrichment = raw$FoldEnrichment,
    Gene_count = raw$Count,
    P_value = raw$pvalue,
    FDR = raw$p.adjust,
    q_value = raw$qvalue,
    Genes = clean_symbol_list(raw$geneID),
    stringsAsFactors = FALSE
  )
}

sample_meta <- s1$sample_meta
sample_meta <- sample_meta[match(colnames(s1$batch_corrected_logcpm), sample_meta$sample_id), , drop = FALSE]
sample_labels <- paste(sample_meta$pcs, sample_meta$condition_label, sep = " | ")
names(sample_labels) <- sample_meta$sample_id

outputs <- list()

for (method_key in names(method_specs)) {
  spec <- method_specs[[method_key]]
  message("Building ", spec$title)
  method_data <- s5$all_method_results[
    as.character(s5$all_method_results$method) == spec$method,
    ,
    drop = FALSE
  ]
  if (nrow(method_data) == 0) stop("No cached results found for ", spec$method)
  method_data$contrast <- as.character(method_data$contrast)
  method_data$contrast_label <- as.character(method_data$contrast_label)
  method_data <- fill_missing_symbols(method_data)
  method_data$has_symbol <- !is.na(method_data$SYMBOL) & trimws(method_data$SYMBOL) != ""
  method_data$is_DEG <- !is.na(method_data$FDR) & method_data$FDR < fdr_threshold & abs(method_data$logFC) >= abs_log2fc_threshold
  method_data$Direction <- ifelse(method_data$is_DEG & method_data$logFC > 0, "Up", ifelse(method_data$is_DEG, "Down", "Not significant"))

  summary <- make_summary(method_data)
  gene_tables <- make_symbol_gene_tables(method_data)
  go_path <- file.path(table_dir, paste0(spec$prefix, "_GO_ORA_all_results.csv"))
  if (!file.exists(go_path)) stop("Required GO result is missing: ", go_path)
  go_terms <- make_go_table(go_path)

  gene_tables$expression$Sample <- unname(sample_labels[gene_tables$expression$Sample])
  sample_order <- unname(sample_labels[colnames(s1$batch_corrected_logcpm)])

  if (any(grepl("ENSG[0-9]+", gene_tables$degs$Symbol, ignore.case = TRUE))) {
    stop("An Ensembl identifier remained in the symbol-only DEG table for ", spec$method)
  }
  if (any(grepl("ENSG[0-9]+", go_terms$Genes, ignore.case = TRUE))) {
    stop("An Ensembl identifier remained in the symbol-only GO gene lists for ", spec$method)
  }
  if (any(duplicated(paste(gene_tables$degs$Contrast_key, gene_tables$degs$Symbol, sep = "||")))) {
    stop("Duplicate symbols remained within a contrast for ", spec$method)
  }

  app_dir <- app_dirs[[method_key]]
  dir.create(app_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(app_template, file.path(app_dir, "app.R"), overwrite = TRUE)

  technical_summary <- paste0(
    spec$method_label, " identified ",
    paste0(summary$Called_DEGs, " DEGs for ", summary$Contrast, collapse = "; "),
    ". Symbol-only tables show ",
    paste0(summary$Symbol_mapped_DEGs, " mapped genes", collapse = " and "),
    "; unmapped records remain included in the called DEG totals but are not displayed as gene labels."
  )

  report_data <- list(
    title = spec$title,
    subtitle = spec$subtitle,
    technical_summary = technical_summary,
    file_prefix = paste0(spec$prefix, "_symbol_only_v2"),
    generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    contrast_order = names(contrast_specs),
    contrast_labels = vapply(contrast_specs, `[[`, character(1), "label"),
    summary = summary,
    genes = gene_tables$degs,
    all_genes = gene_tables$all_genes,
    expression = gene_tables$expression,
    sample_order = sample_order,
    go_terms = go_terms,
    scope_notes = c(
      "Comparisons are 2DG + OX versus control and LDHA + OX versus control.",
      "DEGs require FDR < 0.05 and an absolute fold change of at least 1.5, equivalent to |log2FC| >= log2(1.5).",
      "Positive log2FC indicates higher expression in treatment relative to control.",
      "All user-facing gene labels and gene-list columns use symbols only. Records without a current symbol are excluded from gene displays and downloads."
    ),
    method_notes = c(
      spec$method_note,
      "ComBat-seq used PCS1 through PCS5 as separate batches while preserving treatment condition.",
      "Gene Ontology BP, MF, and CC over-representation analyses use the retained expressed genes with Entrez mappings as the universe and Benjamini-Hochberg adjustment.",
      "Within a contrast, multiple records mapping to the same symbol are resolved by retaining the row with the smallest FDR. No significant duplicate-symbol rows remained in these data."
    ),
    limitation_notes = c(
      "Corrected-count tests do not directly model the paired PCS block after correction and can be more permissive than the conservative edgeR quasi-likelihood model with PCS in the design.",
      "Symbol-only presentation omits significant records that lack a current symbol; called DEG totals therefore exceed the number of rows shown for some contrasts.",
      "GO results depend on annotation coverage, the DEG threshold, the selected DEG direction, and redundancy among ontology terms.",
      "Upstream FASTQ or BAM quality metrics, RNA integrity, alignment rates, strandedness, and library-preparation covariates were not available."
    ),
    next_steps = c(
      "Use the selected-row plots to define manuscript figures only after checking effect direction and corrected expression across all PCS blocks.",
      "Prioritize genes and GO terms that agree between the edgeR exact and DESeq2 ComBat-seq reports and remain biologically coherent with the conservative PCS-adjusted model.",
      "Report the symbol-mapping omission counts alongside DEG totals when these tables are used in supplementary material."
    ),
    provenance_note = paste0(
      "Generated from cached stage 01 QC and ComBat-seq outputs, cached stage 05 method results, and stage 06 GO results on ",
      format(Sys.Date(), "%d %B %Y"), "."
    )
  )
  saveRDS(report_data, file.path(app_dir, "report_data.rds"), compress = "xz")

  launcher <- c(
    "file_arg <- grep('^--file=', commandArgs(FALSE), value = TRUE)",
    "app_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub('^--file=', '', file_arg[[1]]), winslash = '/', mustWork = TRUE)) else normalizePath(getwd(), winslash = '/', mustWork = TRUE)",
    "port <- suppressWarnings(as.integer(Sys.getenv('RNASEQ_REPORT_PORT', '')))",
    "if (is.na(port)) port <- NULL",
    "shiny::runApp(app_dir, host = '127.0.0.1', port = port, launch.browser = TRUE)"
  )
  writeLines(launcher, file.path(app_dir, "start_report.R"))

  manifest <- c(
    paste("title:", spec$title),
    paste("method:", spec$method_label),
    paste("generated:", report_data$generated),
    "threshold: FDR < 0.05 and fold change >= 1.5",
    "batch correction: ComBat-seq with each PCS block treated as a separate batch and condition preserved",
    "gene display policy: symbols only; unmapped significant records omitted from display and reported in summary",
    paste("source DEG cache:", stage5_rds),
    paste("source GO table:", go_path)
  )
  writeLines(manifest, file.path(app_dir, "report_manifest.txt"))

  deg_export <- gene_tables$degs[, c(
    "Symbol", "Gene_name", "Contrast", "Direction", "log2FC", "Fold_change",
    "P_value", "FDR", "Average_expression"
  ), drop = FALSE]
  go_export <- go_terms[, c(
    "GO_ID", "Term", "Ontology", "DEG_set", "Contrast", "Gene_ratio",
    "Background_ratio", "Fold_enrichment", "Gene_count", "P_value", "FDR", "q_value", "Genes"
  ), drop = FALSE]
  summary_export <- summary[, c(
    "Contrast", "Genes_tested", "Called_DEGs", "Symbol_mapped_DEGs", "Unmapped_DEGs",
    "Mapped_Up", "Mapped_Down", "Top_symbol"
  ), drop = FALSE]

  deg_csv <- file.path(table_dir, paste0(spec$prefix, "_symbol_only_DEGs_v2.csv"))
  go_csv <- file.path(table_dir, paste0(spec$prefix, "_symbol_only_GO_ORA_v2.csv"))
  workbook <- file.path(table_dir, paste0(spec$prefix, "_symbol_only_interactive_v2.xlsx"))
  write.csv(deg_export, deg_csv, row.names = FALSE, na = "")
  write.csv(go_export, go_csv, row.names = FALSE, na = "")
  write.xlsx(
    list(
      DEG_summary = summary_export,
      Symbol_DEGs = deg_export,
      Symbol_GO_ORA = go_export,
      Methods = data.frame(
        Item = c("Method", "Batch correction", "DEG threshold", "Symbol policy", "GO universe"),
        Definition = c(
          spec$method_label,
          "ComBat-seq with PCS1 through PCS5 as separate batches and condition preserved",
          "FDR < 0.05 and absolute fold change >= 1.5",
          "Symbols only; unmapped records omitted from gene tables and counted in DEG_summary",
          "Retained expressed genes with Entrez mappings"
        )
      )
    ),
    workbook,
    asTable = TRUE,
    overwrite = TRUE
  )

  outputs[[method_key]] <- list(
    app_dir = app_dir,
    app_file = file.path(app_dir, "app.R"),
    start_file = file.path(app_dir, "start_report.R"),
    workbook = workbook,
    deg_csv = deg_csv,
    go_csv = go_csv,
    summary = summary
  )
}

stage_data <- list(
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  version_tag = version_tag,
  threshold = list(FDR = fdr_threshold, fold_change = fc_threshold, abs_log2FC = abs_log2fc_threshold),
  outputs = outputs
)
saveRDS(stage_data, stage_rds)
write_stage_stamp(stage_name)
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_07_symbol_interactive_reports.txt"))
message("Stage 7 complete: ", stage_rds)
