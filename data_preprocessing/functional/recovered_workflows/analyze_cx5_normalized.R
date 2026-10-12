# Recovered numerical analysis from SSc_fibroblast_analysis_NEW.Rmd.
# Starts at normalized donor/experiment-level percentages, not raw image features.
analyze_cx5_normalized <- function(workbook) {
  wide <- as.data.frame(readxl::read_excel(workbook, sheet = "Normalized_Data"))
  if (!all(c("Channel", "Condition") %in% names(wide))) stop("Channel and Condition required.")
  experiments <- grep("^Exp", names(wide), value = TRUE)
  if (!length(experiments) || !all(vapply(wide[experiments], is.numeric, logical(1))))
    stop("Numeric Exp columns required.")
  channels <- c("aSMA", "Stress fiber")
  conditions <- c("Control", "2DG", "LDHA", "OXPHOS", "OX+2DG", "OX+LDHA")
  wide$Channel[wide$Channel == "StressFiber"] <- "Stress fiber"
  if (any(!wide$Channel %in% channels) || any(!wide$Condition %in% conditions))
    stop("Unexpected channel or condition.")
  if (anyDuplicated(wide[c("Channel", "Condition")])) stop("Duplicate channel/condition rows.")
  summary <- list()
  for (ch in channels) for (condition in conditions) {
    row <- wide[wide$Channel == ch & wide$Condition == condition, experiments, drop = FALSE]
    if (nrow(row) != 1L) stop("Exactly one row per channel/condition required.")
    v <- unlist(row, use.names = FALSE)
    v <- v[!is.na(v)]
    if (!length(v) || any(!is.finite(v))) stop("No finite observations.")
    summary[[length(summary) + 1L]] <- data.frame(Channel = ch, Condition = condition,
      n = length(v), mean = mean(v), sd = sd(v), sem = sd(v) / sqrt(length(v)),
      p.value = if (condition != "Control" && length(v) >= 2L) t.test(v, mu = 100)$p.value else NA_real_)
  }
  summary <- do.call(rbind, summary)
  stats <- summary[summary$Condition != "Control", , drop = FALSE]
  stats$p.adj <- ave(stats$p.value, stats$Channel, FUN = function(p) p.adjust(p, method = "BH"))
  stats$significance <- ifelse(is.na(stats$p.adj), "", ifelse(stats$p.adj < .001, "***",
    ifelse(stats$p.adj < .01, "**", ifelse(stats$p.adj < .05, "*", "ns"))))
  list(summary = summary[setdiff(names(summary), "p.value")], statistics = stats)
}
if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  if (length(args) != 2L) stop("Usage: Rscript analyze_cx5_normalized.R INPUT.xlsx OUTPUT_DIR")
  result <- analyze_cx5_normalized(args[1L])
  dir.create(args[2L], recursive = TRUE, showWarnings = FALSE)
  write.csv(result$summary, file.path(args[2L], "SSc_summary_table.csv"), row.names = FALSE)
  write.csv(result$statistics, file.path(args[2L], "SSc_statistics_table.csv"), row.names = FALSE)
}
