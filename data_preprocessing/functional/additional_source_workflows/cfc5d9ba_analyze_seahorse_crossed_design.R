source("data_preprocessing/common/config.R")
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lmerTest)
  library(emmeans)
})

input_file <- project_path("Revision/Seahorse Ec-fib values.xlsx")
output_dir <- file.path(getwd(), "seahorse_crossed_design_analysis")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

raw <- read_excel(input_file, sheet = "Sheet1", col_names = FALSE, .name_repair = "minimal")

# Preserve the original EC-line order without publishing source donor codes.
ec_source_codes <- strsplit(Sys.getenv("METABOLIC_EC_SOURCE_CODES",
  unset = "EC1,EC2,EC3"), ",", fixed = TRUE)[[1L]]
stopifnot(length(ec_source_codes) == 3L, !anyDuplicated(ec_source_codes))
observed_ec_codes <- unique(as.character(na.omit(raw[[1]][-1])))
if (!all(observed_ec_codes %in% ec_source_codes)) {
  stop("Set METABOLIC_EC_SOURCE_CODES to match the input workbook in EC-line 1, 2, 3 order.")
}

dat <- tibble(
  ec_line_source = raw[[1]][-1],
  fibroblast_observed = as.character(raw[[2]][-1]),
  baseline_ocr = as.numeric(raw[[3]][-1]),
  atp_linked = as.numeric(raw[[5]][-1]),
  proton_leak = as.numeric(raw[[7]][-1]),
  maximal_ocr = as.numeric(raw[[9]][-1]),
  spare_capacity = as.numeric(raw[[11]][-1]),
  non_mito = as.numeric(raw[[13]][-1])
) %>%
  fill(ec_line_source) %>%
  mutate(
    ec_line = factor(
      as.character(ec_line_source),
      levels = ec_source_codes,
      labels = c("1", "2", "3")
    ),
    fibroblast_observed = recode(
      fibroblast_observed,
      "PT.1" = "PT1", "Pt.2" = "PT2", "Pt.3" = "PT3"
    ),
    condition = if_else(grepl("^NH", fibroblast_observed), "Healthy", "Progressive SSc"),
    condition = factor(condition, levels = c("Healthy", "Progressive SSc")),
    # Retain the two healthy fibroblast donors with complementary EC-line coverage as separate healthy fibroblast donors.
    fibroblast_model = factor(fibroblast_observed)
  )

write.csv(dat, file.path(output_dir, "seahorse_analysis_dataset.csv"), row.names = FALSE)

design_counts <- dat %>%
  count(ec_line_source, ec_line, condition, fibroblast_observed,
        fibroblast_model, name = "n")
write.csv(design_counts, file.path(output_dir, "design_counts.csv"), row.names = FALSE)

exact_blocked_donor_permutation <- function(d) {
  donor_condition <- d %>% distinct(fibroblast_model, condition)
  donors <- as.character(donor_condition$fibroblast_model)
  n_healthy <- sum(donor_condition$condition == "Healthy")
  stopifnot(length(donors) == 10L, n_healthy == 4L)

  observed <- unname(coef(lm(value ~ condition + ec_line, data = d))[
    "conditionProgressive SSc"
  ])
  assignments <- combn(donors, n_healthy, simplify = FALSE)
  perm_diffs <- vapply(assignments, function(healthy_donors) {
    d_perm <- d %>%
      mutate(
        perm_condition = factor(
          if_else(as.character(fibroblast_model) %in% healthy_donors,
                  "Healthy", "Progressive SSc"),
          levels = c("Healthy", "Progressive SSc")
        )
      )
    unname(coef(lm(value ~ perm_condition + ec_line, data = d_perm))[
      "perm_conditionProgressive SSc"
    ])
  }, numeric(1))
  p <- mean(abs(perm_diffs) >= abs(observed) - sqrt(.Machine$double.eps))
  c(estimate = observed, p_value = p, permutations = length(perm_diffs))
}

analyze_endpoint <- function(endpoint, label) {
  d <- dat %>%
    select(ec_line, fibroblast_observed, fibroblast_model, condition,
           value = all_of(endpoint))

  full_interaction <- lmer(
    value ~ condition * ec_line + (1 | fibroblast_model),
    data = d,
    REML = TRUE
  )
  full_anova <- anova(full_interaction, ddf = "Kenward-Roger")
  primary <- lmer(
    value ~ condition + ec_line + (1 | fibroblast_model),
    data = d,
    REML = TRUE
  )
  primary_anova <- anova(primary, ddf = "Kenward-Roger")
  primary_emm <- emmeans(primary, ~ condition, lmer.df = "kenward-roger")
  primary_contrast <- as.data.frame(summary(
    contrast(primary_emm, list("Progressive SSc - Healthy" = c(-1, 1))),
    infer = c(TRUE, TRUE), adjust = "none"
  ))
  primary_emm_df <- as.data.frame(primary_emm)

  donor <- d %>%
    group_by(fibroblast_model, condition) %>%
    summarise(value = mean(value), .groups = "drop")
  perm <- exact_blocked_donor_permutation(d)
  wilcox <- wilcox.test(value ~ condition, data = donor, exact = TRUE, conf.int = FALSE)

  naive_welch <- t.test(value ~ condition, data = d)
  naive_pooled <- t.test(value ~ condition, data = d, var.equal = TRUE)
  naive_wilcox <- wilcox.test(value ~ condition, data = d, exact = FALSE)
  blocked_lm <- summary(lm(value ~ condition + ec_line, data = d))$coefficients[
    "conditionProgressive SSc", , drop = TRUE
  ]

  log_row <- tibble(
    log_effect_ratio = NA_real_,
    log_ci_low_ratio = NA_real_,
    log_ci_high_ratio = NA_real_,
    log_lmm_p = NA_real_
  )
  if (all(d$value > 0)) {
    log_fit <- lmer(
      log(value) ~ condition + ec_line + (1 | fibroblast_model),
      data = d,
      REML = TRUE
    )
    log_emm <- emmeans(log_fit, ~ condition, lmer.df = "kenward-roger")
    lc <- as.data.frame(summary(
      contrast(log_emm, list("Progressive SSc - Healthy" = c(-1, 1))),
      infer = c(TRUE, TRUE), adjust = "none"
    ))
    log_row <- tibble(
      log_effect_ratio = exp(lc$estimate),
      log_ci_low_ratio = exp(lc$lower.CL),
      log_ci_high_ratio = exp(lc$upper.CL),
      log_lmm_p = lc$p.value
    )
  }

  residual_shapiro <- shapiro.test(residuals(primary))
  group_summary <- d %>%
    group_by(condition) %>%
    summarise(
      combination_n = n(),
      fibroblast_donors = n_distinct(fibroblast_model),
      raw_mean = mean(value),
      raw_sd = sd(value),
      .groups = "drop"
    )

  result <- tibble(
    endpoint = endpoint,
    label = label,
    estimate_ssc_minus_healthy = primary_contrast$estimate,
    ci_low = primary_contrast$lower.CL,
    ci_high = primary_contrast$upper.CL,
    standard_error = primary_contrast$SE,
    df_kenward_roger = primary_contrast$df,
    lmm_kenward_roger_p = primary_contrast$p.value,
    condition_by_ec_interaction_p = full_anova["condition:ec_line", "Pr(>F)"],
    healthy_adjusted_mean = primary_emm_df$emmean[primary_emm_df$condition == "Healthy"],
    ssc_adjusted_mean = primary_emm_df$emmean[primary_emm_df$condition == "Progressive SSc"],
    random_donor_variance = as.data.frame(VarCorr(primary))$vcov[1],
    singular_fit = lme4::isSingular(primary),
    donor_exact_permutation_p = unname(perm[["p_value"]]),
    donor_exact_wilcoxon_p = wilcox$p.value,
    naive_welch_p = naive_welch$p.value,
    naive_pooled_t_p = naive_pooled$p.value,
    naive_combination_wilcoxon_p = naive_wilcox$p.value,
    ec_blocked_lm_ignoring_donor_p = blocked_lm[["Pr(>|t|)"]],
    residual_shapiro_p = residual_shapiro$p.value
  ) %>% bind_cols(log_row)

  list(result = result, group_summary = group_summary, data = d,
       model = primary, donor = donor)
}

endpoints <- c(
  baseline_ocr = "Basal OCR",
  atp_linked = "ATP-linked respiration",
  proton_leak = "Proton leak",
  maximal_ocr = "Maximal OCR",
  spare_capacity = "Spare respiratory capacity",
  non_mito = "Non-mitochondrial respiration"
)

analyses <- Map(analyze_endpoint, names(endpoints), unname(endpoints))
names(analyses) <- names(endpoints)
results <- bind_rows(lapply(analyses, `[[`, "result"))
results <- results %>%
  mutate(
    holm_p_three_displayed = if_else(
      endpoint %in% c("baseline_ocr", "atp_linked", "maximal_ocr"),
      p.adjust(
        lmm_kenward_roger_p[endpoint %in% c("baseline_ocr", "atp_linked", "maximal_ocr")],
        method = "holm"
      )[match(
        endpoint,
        endpoint[endpoint %in% c("baseline_ocr", "atp_linked", "maximal_ocr")]
      )],
      NA_real_
    ),
    holm_exact_permutation_p_three_displayed = if_else(
      endpoint %in% c("baseline_ocr", "atp_linked", "maximal_ocr"),
      p.adjust(
        donor_exact_permutation_p[endpoint %in% c("baseline_ocr", "atp_linked", "maximal_ocr")],
        method = "holm"
      )[match(
        endpoint,
        endpoint[endpoint %in% c("baseline_ocr", "atp_linked", "maximal_ocr")]
      )],
      NA_real_
    )
  )
group_summaries <- bind_rows(Map(
  function(x, nm) mutate(x$group_summary, endpoint = nm, .before = 1),
  analyses, names(analyses)
))
donor_means <- bind_rows(Map(
  function(x, nm) mutate(x$donor, endpoint = nm, .before = 1),
  analyses, names(analyses)
))

write.csv(results, file.path(output_dir, "statistical_results_all_endpoints.csv"), row.names = FALSE)
write.csv(group_summaries, file.path(output_dir, "group_summaries_all_endpoints.csv"), row.names = FALSE)
write.csv(donor_means, file.path(output_dir, "donor_means_all_endpoints.csv"), row.names = FALSE)

capture.output(
  lapply(analyses, function(x) summary(x$model)),
  file = file.path(output_dir, "primary_lmm_model_summaries.txt")
)

figure_endpoints <- c("baseline_ocr", "atp_linked", "maximal_ocr")
plot_dat <- dat %>%
  select(ec_line, fibroblast_observed, fibroblast_model, condition,
         all_of(figure_endpoints)) %>%
  pivot_longer(all_of(figure_endpoints), names_to = "endpoint", values_to = "value") %>%
  mutate(
    endpoint = factor(endpoint, levels = figure_endpoints,
                      labels = endpoints[figure_endpoints])
  )
plot_stats <- results %>%
  filter(endpoint %in% figure_endpoints) %>%
  mutate(
    endpoint = factor(endpoint, levels = figure_endpoints,
                      labels = endpoints[figure_endpoints]),
    annotation = sprintf("P = %.4f", lmm_kenward_roger_p)
  ) %>%
  left_join(
    plot_dat %>%
      group_by(endpoint) %>%
      summarise(
        annotation_y = max(value) + 0.13 * diff(range(value)),
        .groups = "drop"
      ),
    by = "endpoint"
  )

p <- ggplot(plot_dat, aes(condition, value, color = condition)) +
  geom_boxplot(aes(group = condition), width = 0.52, outlier.shape = NA,
               color = "grey30", fill = NA, linewidth = 0.45) +
  geom_point(aes(shape = ec_line),
             position = position_jitter(width = 0.12, height = 0), size = 2.0,
             alpha = 0.85) +
  facet_wrap(~ endpoint, scales = "free_y", nrow = 1) +
  geom_text(
    data = plot_stats,
    aes(x = 1.5, y = annotation_y, label = annotation),
    inherit.aes = FALSE, vjust = 0, size = 3.2
  ) +
  scale_color_manual(values = c("Healthy" = "#2878B5", "Progressive SSc" = "#E83F3F")) +
  scale_shape_manual(values = c("1" = 16, "2" = 17, "3" = 15)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = NULL, y = "OCR-derived value", color = NULL, shape = "EC line") +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "right",
    panel.spacing = unit(1.0, "lines")
  )

ggsave(file.path(output_dir, "seahorse_primary_endpoints_mixed_model.svg"), p,
       width = 8.4, height = 3.25, units = "in", bg = "white")
ggsave(file.path(output_dir, "seahorse_primary_endpoints_mixed_model.pdf"), p,
       width = 8.4, height = 3.25, units = "in", bg = "white")
ggsave(file.path(output_dir, "seahorse_primary_endpoints_mixed_model.png"), p,
       width = 8.4, height = 3.25, units = "in", dpi = 300, bg = "white")

writeLines(c(
  "Interaction check: outcome ~ condition * EC line + (1 | fibroblast donor).",
  "Because no displayed endpoint had evidence of condition-by-EC interaction,",
  "the primary marginal model is outcome ~ condition + EC line + (1 | fibroblast donor).",
  "P values: two-sided Kenward-Roger tests for the condition marginal contrast.",
  "Sensitivity: exact donor-label permutation test preserving each donor's observed EC-line measurements and adjusting for EC line.",
  "the two healthy fibroblast donors with complementary EC-line coverage are retained as separate healthy fibroblast donors in all statistical tests.",
  "EC source codes are ordered by METABOLIC_EC_SOURCE_CODES and displayed as EC lines 1, 2 and 3, respectively."
), file.path(output_dir, "README.txt"))

print(results, width = Inf)
