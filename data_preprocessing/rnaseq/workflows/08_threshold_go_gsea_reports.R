# Scientific workflow derived from Revision/RNAseq PCS met inhibitors/results/scripts/08_threshold_go_gsea_reports.R
# Run from the repository root with METABOLIC_INPUT_DIR set to the external project data.
source("data_preprocessing/common/config.R")
script_dir <- {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)) else getwd()
}
source(file.path(script_dir, "00_config.R"))

stage_name <- "08_threshold_go_gsea_reports"
stage_rds <- file.path(cache_dir, paste0(stage_name, ".rds"))
stage1_rds <- file.path(cache_dir, "01_qc_batch_correction.rds")
stage5_rds <- file.path(cache_dir, "05_method_threshold_sensitivity.rds")
stage6_rds <- file.path(cache_dir, "06_combatseq_deg_go_reports.rds")
app_template <- file.path(script_dir, "threshold_go_gsea_report_app.R")
version_tag <- "threshold_GO_GSEA_v5"
gsea_signature <- "signed_log10p_entrez_GO_ALL_fgsea_v1"

method_specs <- list(
  edgeR_exact_ComBatSeq = list(
    method = "edgeR_exact_ComBatSeq",
    title = "edgeR exact ComBat-seq threshold, GO, and GSEA report",
    prefix = "edgeR_exact_ComBatSeq_FDR05_FC15",
    method_label = "edgeR exact test on ComBat-seq corrected counts",
    method_note = "edgeR exactTest was run on condition-preserving ComBat-seq PCS-corrected counts as a two-group comparison for each treatment versus control."
  ),
  DESeq2_ComBatSeq = list(
    method = "DESeq2_ComBatSeq",
    title = "DESeq2 ComBat-seq threshold, GO, and GSEA report",
    prefix = "DESeq2_ComBatSeq_FDR05_FC15",
    method_label = "DESeq2 Wald test on ComBat-seq corrected counts",
    method_note = "DESeq2 Wald tests were run on condition-preserving ComBat-seq PCS-corrected counts with design ~ condition for each treatment versus control."
  )
)

app_dirs <- vapply(method_specs, function(spec) {
  file.path(report_dir, paste0(spec$prefix, "_", version_tag))
}, character(1))
expected_outputs <- unlist(lapply(app_dirs, function(path) {
  file.path(path, c("app.R", "report_data.rds", "start_report.R", "report_manifest.txt", "chart_map.csv"))
}))

deps <- c(
  stage1_rds,
  stage5_rds,
  stage6_rds,
  app_template,
  file.path(script_dir, "00_config.R"),
  file.path(script_dir, "08_threshold_go_gsea_reports.R")
)

if (stage_is_fresh(stage_rds, deps) && all(file.exists(expected_outputs))) {
  message("Using cached ", stage_name, ": ", stage_rds)
  quit(save = "no", status = 0)
}

message("Stage 8: threshold-responsive DEG/GO reports with cached GSEA")
s1 <- read_stage(stage1_rds)
s5 <- read_stage(stage5_rds)

gsea_cache_dir <- file.path(cache_dir, "08_method_gsea")
dir.create(gsea_cache_dir, recursive = TRUE, showWarnings = FALSE)

fill_missing_symbols <- function(data) {
  missing <- is.na(data$SYMBOL) | trimws(data$SYMBOL) == ""
  if (!any(missing)) return(data)
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
  if (!is.null(mapped)) data$SYMBOL[missing] <- unname(mapped[data$ensembl_id[missing]])
  data
}

clean_gene_name <- function(x) {
  x[is.na(x)] <- ""
  trimws(x)
}

clean_symbol_tokens <- function(x) {
  all_tokens <- unique(unlist(strsplit(as.character(x), "/", fixed = TRUE)))
  all_tokens <- trimws(all_tokens)
  all_tokens <- all_tokens[nzchar(all_tokens)]
  numeric_tokens <- grepl("^[0-9]+$", all_tokens)
  if (any(numeric_tokens)) {
    mapped <- suppressMessages(
      AnnotationDbi::mapIds(
        org.Hs.eg.db,
        keys = all_tokens[numeric_tokens],
        column = "SYMBOL",
        keytype = "ENTREZID",
        multiVals = "first"
      )
    )
    all_tokens[numeric_tokens] <- unname(mapped[all_tokens[numeric_tokens]])
  }
  all_tokens <- unique(all_tokens[!is.na(all_tokens) & nzchar(all_tokens)])
  all_tokens <- all_tokens[!grepl("^ENSG[0-9]+", all_tokens, ignore.case = TRUE)]
  paste(all_tokens, collapse = ", ")
}

build_ranked_list <- function(method_data, contrast) {
  data <- method_data[
    method_data$contrast == contrast &
      !is.na(method_data$ENTREZID) & nzchar(method_data$ENTREZID) &
      is.finite(method_data$PValue) & is.finite(method_data$logFC),
    ,
    drop = FALSE
  ]
  data$rank_raw <- sign(data$logFC) * -log10(pmax(data$PValue, .Machine$double.xmin))
  data <- data[order(-abs(data$rank_raw), data$PValue, data$ENTREZID), , drop = FALSE]
  data <- data[!duplicated(data$ENTREZID), , drop = FALSE]
  tie_break <- match(data$ENTREZID, sort(unique(data$ENTREZID))) / max(1, nrow(data))
  sign_nonzero <- ifelse(data$rank_raw < 0, -1, 1)
  data$rank_score <- data$rank_raw + sign_nonzero * tie_break * 1e-12
  ranked <- data$rank_score
  names(ranked) <- data$ENTREZID
  sort(ranked, decreasing = TRUE)
}

run_or_load_gsea <- function(method_data, spec, contrast) {
  cache_path <- file.path(gsea_cache_dir, paste0(spec$method, "_", contrast, "_GO_ALL_gsea.rds"))
  source_mtime <- as.numeric(file.info(stage5_rds)$mtime)
  if (!force_run && file.exists(cache_path)) {
    cached <- tryCatch(readRDS(cache_path), error = function(e) NULL)
    if (!is.null(cached) && identical(cached$signature, gsea_signature) &&
        isTRUE(cached$source_mtime >= source_mtime)) {
      message("Using cached GSEA: ", cache_path)
      return(cached$gsea)
    }
  }

  ranked <- build_ranked_list(method_data, contrast)
  if (length(ranked) < 100) stop("Too few ranked Entrez genes for GSEA: ", spec$method, " ", contrast)
  message("Running GO ALL GSEA for ", spec$method, " ", contrast, " (", length(ranked), " ranked genes)")
  set.seed(42)
  gsea <- suppressMessages(
    clusterProfiler::gseGO(
      geneList = ranked,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "ALL",
      exponent = 1,
      minGSSize = 10,
      maxGSSize = 500,
      eps = 0,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      verbose = FALSE,
      seed = TRUE,
      by = "fgsea"
    )
  )
  if (nrow(as.data.frame(gsea)) == 0) stop("GSEA returned no terms for ", spec$method, " ", contrast)
  gsea <- suppressMessages(DOSE::setReadable(gsea, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"))
  saveRDS(
    list(signature = gsea_signature, source_mtime = source_mtime, gsea = gsea),
    cache_path,
    compress = "xz"
  )
  gsea
}

gsea_to_table <- function(gsea, contrast) {
  raw <- as.data.frame(gsea)
  if (nrow(raw) == 0) return(data.frame())
  core_symbols <- vapply(raw$core_enrichment, clean_symbol_tokens, character(1))
  data.frame(
    GO_ID = raw$ID,
    Term = raw$Description,
    Ontology = raw$ONTOLOGY,
    Contrast_key = contrast,
    Contrast = contrast_specs[[contrast]]$label,
    Direction = ifelse(raw$NES >= 0, "Positive", "Negative"),
    Set_size = raw$setSize,
    Enrichment_score = raw$enrichmentScore,
    NES = raw$NES,
    Rank_at_max = raw$rank,
    Leading_edge = raw$leading_edge,
    Core_gene_count = vapply(strsplit(core_symbols, ", ", fixed = TRUE), function(x) sum(nzchar(x)), integer(1)),
    P_value = raw$pvalue,
    FDR = raw$p.adjust,
    q_value = raw$qvalue,
    Core_genes = core_symbols,
    stringsAsFactors = FALSE
  )
}

build_method_records <- function(method_data) {
  method_data$has_symbol <- !is.na(method_data$SYMBOL) & trimws(method_data$SYMBOL) != ""
  analysis_results <- data.frame(
    Contrast_key = method_data$contrast,
    log2FC = method_data$logFC,
    P_value = method_data$PValue,
    FDR = method_data$FDR,
    has_symbol = method_data$has_symbol,
    stringsAsFactors = FALSE
  )

  mapped <- method_data[method_data$has_symbol, , drop = FALSE]
  unique_gene_ids <- unique(mapped$gene_id)
  expression_key <- match(mapped$gene_id, unique_gene_ids)
  gene_results <- data.frame(
    row_key = seq_len(nrow(mapped)),
    expression_key = expression_key,
    Symbol = mapped$SYMBOL,
    Gene_name = clean_gene_name(mapped$GENENAME),
    Entrez_ID = as.character(mapped$ENTREZID),
    Contrast_key = mapped$contrast,
    Contrast = mapped$contrast_label,
    log2FC = mapped$logFC,
    Fold_change = 2^abs(mapped$logFC),
    Statistic = mapped$statistic,
    P_value = mapped$PValue,
    FDR = mapped$FDR,
    Average_expression = mapped$average_expression,
    stringsAsFactors = FALSE
  )

  expression_ids <- unique_gene_ids[unique_gene_ids %in% rownames(s1$batch_corrected_logcpm)]
  expression <- if (length(expression_ids) > 0) {
    do.call(rbind, lapply(seq_along(expression_ids), function(i) {
      gene_id <- expression_ids[[i]]
      data.frame(
        expression_key = match(gene_id, unique_gene_ids),
        Sample_key = colnames(s1$batch_corrected_logcpm),
        logCPM = as.numeric(s1$batch_corrected_logcpm[gene_id, ]),
        stringsAsFactors = FALSE
      )
    }))
  } else {
    data.frame(expression_key = integer(), Sample_key = character(), logCPM = numeric())
  }
  list(analysis_results = analysis_results, gene_results = gene_results, expression = expression)
}

threshold_mask <- function(data, basis = "FDR", cutoff = 0.05, fold_change = 1.5) {
  metric <- if (identical(basis, "P value")) data$P_value else data$FDR
  is.finite(metric) & metric < cutoff & is.finite(data$log2FC) & abs(data$log2FC) >= log2(fold_change)
}

deduplicate_symbols <- function(data, basis = "FDR") {
  metric <- if (identical(basis, "P value")) data$P_value else data$FDR
  data <- data[order(data$Contrast_key, data$Symbol, metric, data$FDR, data$P_value, -abs(data$log2FC)), , drop = FALSE]
  data[!duplicated(paste(data$Contrast_key, data$Symbol, sep = "||")), , drop = FALSE]
}

default_go_table <- function(path) {
  raw <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(raw) == 0) return(data.frame())
  data.frame(
    GO_ID = raw$ID,
    Term = raw$Description,
    Ontology = raw$ont,
    DEG_set = raw$direction,
    Contrast = raw$contrast_label,
    Gene_ratio = raw$GeneRatio,
    Background_ratio = raw$BgRatio,
    Fold_enrichment = raw$FoldEnrichment,
    Gene_count = raw$Count,
    P_value = raw$pvalue,
    FDR = raw$p.adjust,
    q_value = raw$qvalue,
    Genes = gsub("/", ", ", raw$geneID, fixed = TRUE),
    stringsAsFactors = FALSE
  )
}

sample_meta <- s1$sample_meta
sample_meta <- sample_meta[match(colnames(s1$batch_corrected_logcpm), sample_meta$sample_id), , drop = FALSE]
sample_labels <- paste(sample_meta$pcs, sample_meta$condition_label, sep = " | ")
names(sample_labels) <- sample_meta$sample_id

chart_map <- data.frame(
  section = c(
    rep("GO ORA", 4),
    rep("GSEA", 4),
    rep("DEG", 3)
  ),
  analytical_question = c(
    "Which selected DEG-derived terms are enriched and by how much?",
    "Which selected terms have the strongest adjusted evidence?",
    "Which selected terms have the largest fold enrichment?",
    "How many selected DEGs contribute to each term?",
    "Which ranked-list terms are positively or negatively enriched?",
    "How large and in which direction is each selected NES?",
    "How many leading-edge genes drive each selected term?",
    "Where does one selected term accumulate along the ranked gene list after normalization to NES?",
    "What are the selected gene effect sizes?",
    "How do selected genes vary across PCS-corrected samples?",
    "Where do selected genes sit in the all-gene significance landscape?"
  ),
  plot_type = c(
    "Dot plot", "Significance bar", "Fold-enrichment lollipop", "Gene-count bar",
    "NES dot plot", "NES bar", "Leading-edge lollipop", "Enrichment curve",
    "Effect-size bar", "Expression heatmap", "Volcano highlight"
  ),
  fields = c(
    "Term, fold enrichment, gene count, FDR", "Term, FDR, DEG set", "Term, fold enrichment, gene count", "Term, gene count, ontology",
    "Term, NES, set size, FDR", "Term, NES, direction", "Term, core gene count, set size, NES", "Rank, NES-normalized trajectory, hit positions",
    "Symbol, log2FC, direction", "Symbol, sample, corrected logCPM", "Symbol, log2FC, FDR"
  ),
  palette = c(rep("Blue-orange plus neutrals", 11)),
  export_contract = c(rep("User-defined 40-300 mm width and height; adaptive typography and wrapped labels", 11)),
  stringsAsFactors = FALSE
)

outputs <- list()

for (method_key in names(method_specs)) {
  spec <- method_specs[[method_key]]
  message("Preparing ", spec$title)
  method_data <- s5$all_method_results[
    as.character(s5$all_method_results$method) == spec$method,
    ,
    drop = FALSE
  ]
  if (nrow(method_data) == 0) stop("No cached results found for ", spec$method)
  method_data$contrast <- as.character(method_data$contrast)
  method_data$contrast_label <- as.character(method_data$contrast_label)
  method_data <- fill_missing_symbols(method_data)
  records <- build_method_records(method_data)

  gsea_objects <- list()
  gsea_tables <- list()
  for (contrast in names(contrast_specs)) {
    gsea_objects[[contrast]] <- run_or_load_gsea(method_data, spec, contrast)
    gsea_tables[[contrast]] <- gsea_to_table(gsea_objects[[contrast]], contrast)
  }
  gsea_results <- do.call(rbind, gsea_tables)
  rownames(gsea_results) <- NULL

  default_analysis_mask <- threshold_mask(records$analysis_results, "FDR", 0.05, 1.5)
  default_gene_mask <- threshold_mask(records$gene_results, "FDR", 0.05, 1.5)
  default_degs <- deduplicate_symbols(records$gene_results[default_gene_mask, , drop = FALSE], "FDR")
  default_degs$Direction <- ifelse(default_degs$log2FC >= 0, "Up", "Down")
  default_summary <- do.call(rbind, lapply(names(contrast_specs), function(contrast) {
    analysis_rows <- records$analysis_results$Contrast_key == contrast
    gene_rows <- default_degs$Contrast_key == contrast
    data.frame(
      Contrast = contrast_specs[[contrast]]$label,
      Genes_tested = sum(analysis_rows),
      Called_DEGs = sum(default_analysis_mask & analysis_rows),
      Symbol_mapped_DEGs = sum(gene_rows),
      Unmapped_or_duplicate_DEGs = sum(default_analysis_mask & analysis_rows) - sum(gene_rows),
      Mapped_Up = sum(gene_rows & default_degs$Direction == "Up"),
      Mapped_Down = sum(gene_rows & default_degs$Direction == "Down"),
      stringsAsFactors = FALSE
    )
  }))

  default_go_path <- file.path(table_dir, paste0(spec$prefix, "_GO_ORA_all_results.csv"))
  default_go <- default_go_table(default_go_path)
  app_dir <- app_dirs[[method_key]]
  dir.create(app_dir, recursive = TRUE, showWarnings = FALSE)
  app_cache_dir <- file.path(app_dir, "analysis_cache")
  dir.create(app_cache_dir, recursive = TRUE, showWarnings = FALSE)
  previous_cache_dir <- file.path(report_dir, paste0(spec$prefix, "_threshold_GO_GSEA_v4"), "analysis_cache")
  if (dir.exists(previous_cache_dir)) {
    previous_cache_files <- list.files(previous_cache_dir, full.names = TRUE)
    if (length(previous_cache_files) > 0) file.copy(previous_cache_files, app_cache_dir, overwrite = FALSE)
  }
  file.copy(app_template, file.path(app_dir, "app.R"), overwrite = TRUE)
  write.csv(chart_map, file.path(app_dir, "chart_map.csv"), row.names = FALSE)

  report_data <- list(
    title = spec$title,
    subtitle = "Threshold-responsive differential expression and GO ORA with threshold-free GO GSEA for 2DG + OX and LDHA + OX versus control.",
    method_key = method_key,
    method_label = spec$method_label,
    method_note = spec$method_note,
    file_prefix = paste0(spec$prefix, "_", version_tag),
    generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    contrast_order = names(contrast_specs),
    contrast_labels = vapply(contrast_specs, `[[`, character(1), "label"),
    default_threshold = list(basis = "FDR", cutoff = 0.05, fold_change = 1.5),
    default_summary = default_summary,
    analysis_results = records$analysis_results,
    gene_results = records$gene_results,
    expression = records$expression,
    sample_labels = sample_labels,
    sample_order = unname(sample_labels[colnames(s1$batch_corrected_logcpm)]),
    universe_entrez = unique(na.omit(records$gene_results$Entrez_ID[nzchar(records$gene_results$Entrez_ID)])),
    gsea_results = gsea_results,
    gsea_objects = gsea_objects,
    gsea_rank_note = "All mapped tested genes were ranked by sign(log2FC) x -log10(raw P value); GO BP, MF, and CC gene sets were tested with fgsea through clusterProfiler.",
    scope_notes = c(
      "DEG thresholds can use either Benjamini-Hochberg FDR or raw P value, together with a minimum absolute fold change.",
      "GO over-representation analysis is recomputed from the exact active DEG set and cached by method, contrast, significance basis, cutoff, and fold-change threshold.",
      "GSEA is intentionally threshold-free and does not change when the DEG threshold changes.",
      "GSEA enrichment trajectories are rescaled so the curve extremum equals the reported normalized enrichment score (NES).",
      "PDF and 600 dpi PNG export width and height are independently adjustable from 40 to 300 mm for every plot.",
      "All displayed gene labels, GO gene lists, and GSEA leading-edge gene lists use symbols only."
    ),
    limitation_notes = c(
      "Raw P-value thresholds do not control the expected false discovery rate and should be labeled exploratory.",
      "Corrected-count tests do not directly model the paired PCS block after ComBat-seq correction and can be more permissive than a PCS-adjusted quasi-likelihood model.",
      "GO ORA is threshold-dependent, whereas GSEA depends on the complete ranked list and its ranking metric.",
      "GO terms are redundant; related terms should not be interpreted as independent biological findings.",
      "Upstream FASTQ or BAM quality metrics, RNA integrity, alignment rates, strandedness, and library-preparation covariates were not available."
    )
  )
  saveRDS(report_data, file.path(app_dir, "report_data.rds"), compress = "xz")

  launcher <- c(
    "launcher_files <- unique(unlist(lapply(sys.frames(), function(frame) {",
    "  value <- frame$ofile",
    "  if (is.null(value)) character(0) else as.character(value)",
    "  })))",
    "launcher_files <- launcher_files[nzchar(launcher_files) & basename(launcher_files) == 'start_report.R']",
    "candidate_dirs <- unique(dirname(normalizePath(launcher_files, winslash = '/', mustWork = FALSE)))",
    "candidate_dirs <- candidate_dirs[file.exists(file.path(candidate_dirs, 'app.R'))]",
    sprintf(
      "fallback_app_dir <- %s",
      encodeString(normalizePath(app_dir, winslash = "/", mustWork = TRUE), quote = "\"")
    ),
    "app_dir <- if (length(candidate_dirs) > 0) candidate_dirs[[1]] else fallback_app_dir",
    "if (!file.exists(file.path(app_dir, 'app.R'))) stop('Could not find app.R next to start_report.R: ', app_dir)",
    "port <- suppressWarnings(as.integer(Sys.getenv('RNASEQ_REPORT_PORT', '')))",
    "if (is.na(port)) port <- NULL",
    "Sys.setenv(RNASEQ_REPORT_APP_DIR = app_dir)",
    "shiny::runApp(app_dir, host = '127.0.0.1', port = port, launch.browser = TRUE)"
  )
  writeLines(launcher, file.path(app_dir, "start_report.R"))
  writeLines(
    c(
      paste("title:", spec$title),
      paste("method:", spec$method_label),
      paste("generated:", report_data$generated),
      "default DEG threshold: FDR < 0.05 and fold change >= 1.5",
      "GO ORA: recomputed from the active threshold-selected DEG set and cached on disk",
      paste("GSEA:", report_data$gsea_rank_note),
      "figure dimensions: user-defined width and height from 40 to 300 mm for PDF and 600 dpi PNG exports",
      paste("source DEG cache:", stage5_rds),
      paste("source corrected expression cache:", stage1_rds)
    ),
    file.path(app_dir, "report_manifest.txt")
  )

  deg_export <- default_degs[, c(
    "Symbol", "Gene_name", "Contrast", "Direction", "log2FC", "Fold_change",
    "P_value", "FDR", "Average_expression"
  ), drop = FALSE]
  gsea_export <- gsea_results[, c(
    "GO_ID", "Term", "Ontology", "Contrast", "Direction", "Set_size", "NES",
    "P_value", "FDR", "q_value", "Core_gene_count", "Core_genes"
  ), drop = FALSE]
  workbook <- file.path(table_dir, paste0(spec$prefix, "_", version_tag, ".xlsx"))
  write.xlsx(
    list(
      Default_DEG_summary = default_summary,
      Default_symbol_DEGs = deg_export,
      Default_GO_ORA = default_go,
      GO_GSEA = gsea_export,
      Methods = data.frame(
        Item = c("Method", "Batch correction", "Default DEG threshold", "Dynamic GO", "GSEA ranking", "Figure dimensions"),
        Definition = c(
          spec$method_label,
          "ComBat-seq with PCS1 through PCS5 as separate batches and condition preserved",
          "FDR < 0.05 and absolute fold change >= 1.5",
          "GO ORA is recomputed from the active DEG selection and cached",
          report_data$gsea_rank_note,
          "Independent user-defined width and height from 40 to 300 mm for PDF and 600 dpi PNG exports"
        )
      )
    ),
    workbook,
    asTable = TRUE,
    overwrite = TRUE
  )
  gsea_csv <- file.path(table_dir, paste0(spec$prefix, "_GO_GSEA_symbol_only_v5.csv"))
  write.csv(gsea_export, gsea_csv, row.names = FALSE, na = "")

  outputs[[method_key]] <- list(
    app_dir = app_dir,
    app_file = file.path(app_dir, "app.R"),
    start_file = file.path(app_dir, "start_report.R"),
    workbook = workbook,
    gsea_csv = gsea_csv,
    default_summary = default_summary,
    gsea_rows = nrow(gsea_results)
  )
}

stage_data <- list(
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  version_tag = version_tag,
  gsea_signature = gsea_signature,
  outputs = outputs
)
saveRDS(stage_data, stage_rds)
write_stage_stamp(stage_name)
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_08_threshold_go_gsea_reports.txt"))
message("Stage 8 complete: ", stage_rds)
