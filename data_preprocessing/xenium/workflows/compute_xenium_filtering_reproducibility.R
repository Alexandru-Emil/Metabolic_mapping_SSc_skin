# Scientific workflow derived from Revision/Nan annotation/Data/RMD files/reviewer_QC_analyses_20260729/03_Xenium_filtering_reproducibility/compute_xenium_filtering_reproducibility.R
# Run from the repository root with METABOLIC_INPUT_DIR set to the external project data.
source("data_preprocessing/common/config.R")
options(stringsAsFactors = FALSE, width = 220)

suppressPackageStartupMessages({
  library(SpatialFeatureExperiment)
  library(SummarizedExperiment)
  library(dplyr)
})

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || is.na(x[1])) y else x
}

script_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
script_dir <- normalizePath(
  dirname(script_file %||% "xenium_filtering_reproducibility/compute_xenium_filtering_reproducibility.R"),
  winslash = "/",
  mustWork = FALSE
)
output_dir <- file.path(script_dir, "outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

project_root <- paste0(
  project_path("Revision")
)
xenium_object_root <- file.path(
  project_root,
  "Nan annotation/Data/Xenium sfe objects"
)

source_paths <- c(
  original_qc_report = file.path(project_root, "Xenium/sfe_QC.html"),
  six_donor_pre_qc = file.path(
    xenium_object_root,
    "Xenium sfe object/sfe object250723154452.rds"
  ),
  six_donor_post_qc = file.path(
    xenium_object_root,
    "Xenium sfe object/sfe_afterqc.rds"
  ),
  ssc230_pre_join = file.path(
    xenium_object_root,
    "XeniumSfe_ExcludedValidation_260729115149.rds"
  ),
  ssc230_construction_report = file.path(
    xenium_object_root,
    "XeniumSfe_ExcludedValidation_20260729_102018.html"
  )
)
stopifnot(all(file.exists(source_paths)))

# Literal thresholds recovered from the executed 2025-07-23 QC report.
six_thresholds <- data.frame(
  sample = c("Validation5", "Validation6", "Validation1", "Validation2", "Validation3", "Validation4"),
  area_min_um2_inclusive = c(15, 15, 19, 13, 13, 14),
  area_max_um2_inclusive = c(270, 250, 220, 300, 280, 300),
  stringsAsFactors = FALSE
)

pre <- readRDS(source_paths[["six_donor_pre_qc"]])
post <- readRDS(source_paths[["six_donor_post_qc"]])
ssc230 <- readRDS(source_paths[["ssc230_pre_join"]])

pre_cd <- as.data.frame(colData(pre), check.names = FALSE)
post_cd <- as.data.frame(colData(post), check.names = FALSE)
ssc230_cd <- as.data.frame(colData(ssc230), check.names = FALSE)

ssc230_area_min <- 5
ssc230_area_max <- 400
ssc230_area_pass <- with(
  ssc230_cd,
  cell_area >= ssc230_area_min & cell_area <= ssc230_area_max
)
stopifnot(all(ssc230_area_pass))

count_by_sample <- function(metadata) {
  data.frame(sample = as.character(metadata$SampleId)) |>
    count(sample, name = "n_cells")
}

six_counts <- six_thresholds |>
  left_join(
    count_by_sample(pre_cd) |> rename(cells_before_qc = n_cells),
    by = "sample"
  ) |>
  left_join(
    count_by_sample(post_cd) |> rename(cells_retained = n_cells),
    by = "sample"
  ) |>
  mutate(
    cells_excluded = cells_before_qc - cells_retained,
    percent_retained = 100 * cells_retained / cells_before_qc,
    transcript_threshold = "> 20",
    detected_gene_threshold = "> 10",
    negative_control_threshold = "< 0.001",
    area_filter = paste0(
      area_min_um2_inclusive,
      "-",
      area_max_um2_inclusive,
      " um^2 (inclusive)"
    ),
    negative_control_formula = paste0(
      "(control_probe_counts + genomic_control_counts + ",
      "control_codeword_counts + unassigned_codeword_counts) / ",
      "(transcript_counts + the same four control-count terms)"
    ),
    pipeline = "Original six-donor QC"
  )

ssc230_counts <- data.frame(
  sample = "ExcludedValidation",
  area_min_um2_inclusive = ssc230_area_min,
  area_max_um2_inclusive = ssc230_area_max,
  cells_before_qc = nrow(ssc230_cd),
  cells_retained = sum(as.logical(ssc230_cd$qc_pass), na.rm = TRUE),
  stringsAsFactors = FALSE
) |>
  mutate(
    cells_excluded = cells_before_qc - cells_retained,
    percent_retained = 100 * cells_retained / cells_before_qc,
    transcript_threshold = "> 20",
    detected_gene_threshold = "> 10",
    negative_control_threshold = "< 0.001",
    area_filter = paste0(
      ssc230_area_min,
      "-",
      ssc230_area_max,
      " um^2 (inclusive; excludes 0 cells)"
    ),
    negative_control_formula = paste0(
      "(control_probe_counts + control_codeword_counts + ",
      "unassigned_codeword_counts) / total_counts"
    ),
    pipeline = "Later ExcludedValidation QC with explicit non-excluding area bounds"
  )

threshold_retention <- bind_rows(six_counts, ssc230_counts) |>
  select(
    sample,
    pipeline,
    transcript_threshold,
    detected_gene_threshold,
    negative_control_threshold,
    area_filter,
    area_min_um2_inclusive,
    area_max_um2_inclusive,
    cells_before_qc,
    cells_retained,
    cells_excluded,
    percent_retained,
    negative_control_formula
  )

six_control_numerator <- with(
  post_cd,
  control_probe_counts +
    genomic_control_counts +
    control_codeword_counts +
    unassigned_codeword_counts
)
six_negative_recomputed <- six_control_numerator /
  (post_cd$transcript_counts + six_control_numerator)

ssc230_negative_recomputed <- with(
  ssc230_cd,
  (control_probe_counts +
    control_codeword_counts +
    unassigned_codeword_counts) / total_counts
)

ssc230_six_donor_numerator <- with(
  ssc230_cd,
  control_probe_counts +
    genomic_control_counts +
    control_codeword_counts +
    unassigned_codeword_counts
)
ssc230_negative_six_donor_formula <- ssc230_six_donor_numerator /
  (ssc230_cd$transcript_counts + ssc230_six_donor_numerator)
ssc230_pass_current <- as.logical(ssc230_cd$qc_pass)
ssc230_base_pass <- with(
  ssc230_cd,
  transcript_counts > 20 & nFeature > 10
)
ssc230_pass_six_donor_negative_formula <- ssc230_base_pass &
  ssc230_negative_six_donor_formula < 0.001 &
  ssc230_area_pass

ssc230_area_audit <- data.frame(
  sample = "ExcludedValidation",
  observed_area_min_um2 = min(ssc230_cd$cell_area, na.rm = TRUE),
  observed_area_max_um2 = max(ssc230_cd$cell_area, na.rm = TRUE),
  applied_area_min_um2_inclusive = ssc230_area_min,
  applied_area_max_um2_inclusive = ssc230_area_max,
  cells_evaluated = nrow(ssc230_cd),
  cells_passing_area = sum(ssc230_area_pass, na.rm = TRUE),
  cells_excluded_by_area = sum(!ssc230_area_pass, na.rm = TRUE),
  selection_note = paste0(
    "Round outward bounds around the observed area range; ",
    "the area predicate is explicit and reproducible but non-excluding."
  ),
  stringsAsFactors = FALSE
)

post_threshold_check <- post_cd |>
  transmute(
    sample = as.character(SampleId),
    transcript_pass = transcript_counts > 20,
    detected_gene_pass = nFeature > 10,
    negative_control_pass = new_prop_neg_count < 0.001,
    cell_area = cell_area
  ) |>
  left_join(six_thresholds, by = "sample") |>
  mutate(
    area_pass = cell_area >= area_min_um2_inclusive &
      cell_area <= area_max_um2_inclusive
  )

validation <- data.frame(
  check = c(
    "Six-donor post-QC cells are an exact subset of pre-QC cells",
    "No duplicated cell identifiers in six-donor pre-QC object",
    "No duplicated cell identifiers in six-donor post-QC object",
    "All retained six-donor cells have transcript_counts > 20",
    "All retained six-donor cells have nFeature > 10",
    "All retained six-donor cells have negative-control proportion < 0.001",
    "All retained six-donor cells satisfy their inclusive area bounds",
    "Stored six-donor negative-control proportions match the recovered formula",
    "Stored ExcludedValidation negative-control proportions match its construction report",
    "All supplied ExcludedValidation cells satisfy the inclusive 5-400 um^2 bounds",
    "Adding the 5-400 um^2 area predicate leaves ExcludedValidation qc_pass unchanged"
  ),
  passed = c(
    all(colnames(post) %in% colnames(pre)),
    anyDuplicated(colnames(pre)) == 0,
    anyDuplicated(colnames(post)) == 0,
    all(post_threshold_check$transcript_pass),
    all(post_threshold_check$detected_gene_pass),
    all(post_threshold_check$negative_control_pass),
    all(post_threshold_check$area_pass),
    max(
      abs(post_cd$new_prop_neg_count - six_negative_recomputed),
      na.rm = TRUE
    ) < 1e-14,
    isTRUE(all.equal(
      ssc230_cd$new_prop_neg_count,
      ssc230_negative_recomputed,
      tolerance = 1e-15,
      check.attributes = FALSE
    )),
    all(ssc230_area_pass),
    identical(
      ssc230_pass_current,
      with(
        ssc230_cd,
        transcript_counts > 20 &
          nFeature > 10 &
          new_prop_neg_count < 0.001 &
          ssc230_area_pass
      )
    )
  ),
  stringsAsFactors = FALSE
)

formula_sensitivity <- data.frame(
  sample = "ExcludedValidation",
  cells_before_qc = nrow(ssc230_cd),
  retained_current_ssc230_formula = sum(ssc230_pass_current),
  retained_with_original_six_donor_negative_formula_and_5_400_area_filter =
    sum(ssc230_pass_six_donor_negative_formula),
  retained_only_by_current_ssc230_formula =
    sum(ssc230_pass_current & !ssc230_pass_six_donor_negative_formula),
  retained_only_by_original_six_donor_formula =
    sum(!ssc230_pass_current & ssc230_pass_six_donor_negative_formula),
  stringsAsFactors = FALSE
)

blinding_audit <- data.frame(
  audit_item = c(
    "Sample identity available during QC",
    "Clinical progressive/stable phenotype available in pre-QC colData",
    "Downstream cell-phenotype annotations available in pre-QC colData",
    "Formal blinding statement in the executed QC report",
    "Thresholding occurred before ResolVI/Leiden phenotyping"
  ),
  evidence = c(
    paste(
      c("SampleId", "Patient", "Group")[
        c("SampleId", "Patient", "Group") %in% colnames(pre_cd)
      ],
      collapse = ", "
    ),
    as.character(any(grepl(
      "progress|stable|phenotype",
      colnames(pre_cd),
      ignore.case = TRUE
    ))),
    as.character(any(grepl(
      "anno|cluster|leiden|celltype|cell_type",
      colnames(pre_cd),
      ignore.case = TRUE
    ))),
    "No",
    "Yes; the pre-QC object contains no downstream annotation fields"
  ),
  interpretation = c(
    "The analyst was not blinded to sample identity.",
    "Progressive/stable status is not present in the archived pre-QC metadata.",
    "Cell-phenotype labels were unavailable and could not guide threshold selection.",
    "Formal blinding to clinical phenotype should not be claimed.",
    "State that QC was performed before cell phenotyping."
  ),
  stringsAsFactors = FALSE
)

write.csv(
  threshold_retention,
  file.path(output_dir, "xenium_qc_thresholds_and_retention.csv"),
  row.names = FALSE
)
write.csv(
  validation,
  file.path(output_dir, "xenium_qc_validation_checks.csv"),
  row.names = FALSE
)
write.csv(
  formula_sensitivity,
  file.path(output_dir, "ssc230_negative_control_formula_sensitivity.csv"),
  row.names = FALSE
)
write.csv(
  ssc230_area_audit,
  file.path(output_dir, "ssc230_area_threshold_audit.csv"),
  row.names = FALSE
)
write.csv(
  blinding_audit,
  file.path(output_dir, "xenium_qc_blinding_audit.csv"),
  row.names = FALSE
)

source_inventory <- data.frame(
  source = names(source_paths),
  path = unname(source_paths),
  file_size_bytes = unname(file.info(source_paths)$size),
  modified = format(
    unname(file.info(source_paths)$mtime),
    "%Y-%m-%d %H:%M:%S"
  ),
  md5 = unname(tools::md5sum(source_paths)),
  stringsAsFactors = FALSE
)
write.csv(
  source_inventory,
  file.path(output_dir, "source_inventory_md5.csv"),
  row.names = FALSE
)

six_total_before <- sum(six_counts$cells_before_qc)
six_total_retained <- sum(six_counts$cells_retained)
seven_total_before <- six_total_before + ssc230_counts$cells_before_qc
seven_total_retained <- six_total_retained + ssc230_counts$cells_retained

summary_lines <- c(
  paste0(
    "Original six donors: ",
    format(six_total_retained, big.mark = ","),
    " of ",
    format(six_total_before, big.mark = ","),
    " cells retained (",
    sprintf("%.2f", 100 * six_total_retained / six_total_before),
    "%)."
  ),
  paste0(
    "ExcludedValidation: ",
    format(ssc230_counts$cells_retained, big.mark = ","),
    " of ",
    format(ssc230_counts$cells_before_qc, big.mark = ","),
    " cells retained (",
    sprintf(
      "%.2f",
      100 * ssc230_counts$cells_retained / ssc230_counts$cells_before_qc
    ),
    "%). The inclusive 5-400 um^2 area predicate retained all ",
    format(ssc230_area_audit$cells_passing_area, big.mark = ","),
    " supplied ExcludedValidation cells and therefore did not alter the final set."
  ),
  paste0(
    "Current seven-donor joined total: ",
    format(seven_total_retained, big.mark = ","),
    " of ",
    format(seven_total_before, big.mark = ","),
    " supplied cells retained (",
    sprintf("%.2f", 100 * seven_total_retained / seven_total_before),
    "%)."
  ),
  paste0(
    "Applying the original six-donor negative-control formula to ExcludedValidation ",
    "together with the inclusive 5-400 um^2 area predicate would retain ",
    formula_sensitivity$retained_with_original_six_donor_negative_formula_and_5_400_area_filter,
    " cells, ",
    formula_sensitivity$retained_only_by_current_ssc230_formula,
    " fewer than the current ExcludedValidation rule."
  )
)
writeLines(summary_lines, file.path(output_dir, "audit_summary.txt"))

methods_six <- paste0(
  "Cells were retained when they contained more than 20 gene-expression ",
  "transcripts, more than 10 detected genes, a negative-control count ",
  "proportion below 0.001, and a cell area within the following ",
  "sample-specific inclusive bounds: Validation5, 15-270 um^2; Validation6, ",
  "15-250 um^2; Validation1, 19-220 um^2; Validation2, 13-300 um^2; Validation3, ",
  "13-280 um^2; and Validation4, 14-300 um^2. The negative-control ",
  "proportion was calculated as the sum of control-probe, ",
  "genomic-control, control-codeword, and unassigned-codeword counts ",
  "divided by the sum of these counts and gene-expression transcript ",
  "counts. These criteria retained ",
  format(six_total_retained, big.mark = ","),
  " of ",
  format(six_total_before, big.mark = ","),
  " segmented cells: Validation5, ",
  format(
    six_counts$cells_retained[six_counts$sample == "Validation5"],
    big.mark = ","
  ),
  "; Validation6, ",
  format(
    six_counts$cells_retained[six_counts$sample == "Validation6"],
    big.mark = ","
  ),
  "; Validation1, ",
  format(
    six_counts$cells_retained[six_counts$sample == "Validation1"],
    big.mark = ","
  ),
  "; Validation2, ",
  format(
    six_counts$cells_retained[six_counts$sample == "Validation2"],
    big.mark = ","
  ),
  "; Validation3, ",
  format(
    six_counts$cells_retained[six_counts$sample == "Validation3"],
    big.mark = ","
  ),
  "; and Validation4, ",
  format(
    six_counts$cells_retained[six_counts$sample == "Validation4"],
    big.mark = ","
  ),
  ". Thresholds were selected from sample-wise QC distributions before ",
  "ResolVI integration, Leiden clustering, and cell-phenotype assignment. ",
  "Sample identity was visible during QC; however, downstream ",
  "cell-phenotype labels and progressive/stable status were not present ",
  "in the QC object and were not used for threshold selection. ResolVI ",
  "used sample ID as the batch variable. Expression was normalized with ",
  "the Seurat NormalizeData() function using a scale factor of 1,000, ",
  "and cell phenotypes were assigned by Leiden clustering on the ",
  "ResolVI embeddings."
)

methods_seven <- paste0(
  methods_six,
  " ExcludedValidation was processed subsequently using more than 20 ",
  "gene-expression transcripts, more than 10 detected genes, and a ",
  "negative-control proportion below 0.001, together with an inclusive ",
  "cell-area range of 5-400 um^2. The observed ExcludedValidation cell-area range ",
  "was ",
  sprintf("%.3f", ssc230_area_audit$observed_area_min_um2),
  "-",
  sprintf("%.3f", ssc230_area_audit$observed_area_max_um2),
  " um^2, so the area criterion retained all ",
  format(ssc230_area_audit$cells_evaluated, big.mark = ","),
  " supplied ExcludedValidation cells and did not change the final retained set. ",
  "For ExcludedValidation, the negative-control proportion ",
  "was calculated as the sum of negative-control-probe, ",
  "negative-control-codeword, and unassigned-codeword counts divided ",
  "by total counts. The combined filters retained ",
  format(seven_total_retained, big.mark = ","),
  " of ",
  format(seven_total_before, big.mark = ","),
  " supplied cells across seven donors, including ",
  format(ssc230_counts$cells_retained, big.mark = ","),
  " of ",
  format(ssc230_counts$cells_before_qc, big.mark = ","),
  " ExcludedValidation cells."
)

reviewer_response <- paste0(
  "To improve reproducibility, we replaced the phrase 'manually ",
  "adjusted' with the exact sample-specific filtering criteria and now ",
  "report the number of cells retained for each donor. For the original ",
  "six donors, cells were required to have >20 gene-expression ",
  "transcripts, >10 detected genes, a negative-control count proportion ",
  "<0.001, and an area within the following inclusive ranges: Validation5, ",
  "15-270 um^2; Validation6, 15-250 um^2; Validation1, 19-220 um^2; Validation2, ",
  "13-300 um^2; Validation3, 13-280 um^2; and Validation4, 14-300 um^2. These ",
  "criteria retained ",
  format(six_total_retained, big.mark = ","),
  " of ",
  format(six_total_before, big.mark = ","),
  " cells. Threshold selection was performed with sample identity ",
  "visible but before ResolVI integration and Leiden cell phenotyping; ",
  "downstream cell-phenotype labels and progressive/stable status were ",
  "not present in the QC object and were not used to set thresholds. ",
  "For ExcludedValidation, an inclusive area range of 5-400 um^2 was recorded. ",
  "The observed range was ",
  sprintf("%.3f", ssc230_area_audit$observed_area_min_um2),
  "-",
  sprintf("%.3f", ssc230_area_audit$observed_area_max_um2),
  " um^2; consequently, this explicit area predicate retained all ",
  format(ssc230_area_audit$cells_evaluated, big.mark = ","),
  " supplied cells and excluded none. ExcludedValidation retained ",
  format(ssc230_counts$cells_retained, big.mark = ","),
  " cells after the transcript, detected-gene, and negative-control ",
  "criteria. Its negative-control proportion was calculated using the ",
  "subsequent ExcludedValidation construction formula documented in the Methods."
)

writeLines(
  methods_six,
  file.path(output_dir, "methods_text_verified_original_six_donors.txt")
)
writeLines(
  methods_seven,
  file.path(output_dir, "methods_text_exact_current_seven_donors.txt")
)
writeLines(
  reviewer_response,
  file.path(output_dir, "reviewer_response_reproducibility.txt")
)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

if (!all(validation$passed)) {
  stop("At least one QC validation check failed; inspect xenium_qc_validation_checks.csv")
}

message("Wrote Xenium filtering reproducibility outputs to ", output_dir)
