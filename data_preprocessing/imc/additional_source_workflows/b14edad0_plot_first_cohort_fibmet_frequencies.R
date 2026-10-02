source("data_preprocessing/common/config.R")
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

spe_file <- paste0(
  paste0(project_path("Revision/"), "/"),
  "Nan annotation/Data/IMC spe objects/spe_old.rds"
)
output_dir <- file.path(getwd(), "first_cohort_fibmet_frequency_outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

spe <- readRDS(spe_file)
metadata <- as.data.frame(colData(spe))

donor_levels <- paste0("SSc", seq_len(9))
state_levels <- c("Met_hi_Fib", "Other_Fib")
status_levels <- c("progressive", "stable")

state_colors <- c("Met_hi_Fib" = "#E15759", "Other_Fib" = "#4E79A7")
status_colors <- c("progressive" = "#8C2D4F", "stable" = "#1B8A7A")

fibmet_cells <- metadata %>%
  filter(
    Disease == "SSc",
    grepl("^stromal_[1-8]$", as.character(celltype))
  ) %>%
  transmute(
    donor = as.character(Donor),
    patient_id = as.character(patient_id),
    progression_status = as.character(Progression_skin),
    reference_cluster = as.character(celltype),
    metfiblabel = if_else(
      reference_cluster == "stromal_2",
      "Met_hi_Fib",
      "Other_Fib"
    )
  )

unexpected_donors <- setdiff(unique(fibmet_cells$donor), donor_levels)
if (length(unexpected_donors) > 0) {
  stop("Unexpected SSc donor labels: ", paste(unexpected_donors, collapse = ", "))
}

unexpected_status <- setdiff(unique(fibmet_cells$progression_status), status_levels)
if (length(unexpected_status) > 0) {
  stop("Unexpected progression labels: ", paste(unexpected_status, collapse = ", "))
}

donor_annotation <- fibmet_cells %>%
  distinct(donor, patient_id, progression_status) %>%
  mutate(
    donor = factor(donor, levels = donor_levels),
    progression_status = factor(progression_status, levels = status_levels)
  ) %>%
  arrange(donor)

if (any(table(donor_annotation$donor) != 1)) {
  stop("Each donor must have exactly one progression annotation.")
}

frequency_table <- fibmet_cells %>%
  count(donor, patient_id, progression_status, metfiblabel, name = "n_cells") %>%
  group_by(donor) %>%
  mutate(
    n_fibroblasts = sum(n_cells),
    frequency = n_cells / n_fibroblasts
  ) %>%
  ungroup() %>%
  mutate(
    donor = factor(donor, levels = donor_levels),
    metfiblabel = factor(metfiblabel, levels = state_levels),
    progression_status = factor(progression_status, levels = status_levels)
  ) %>%
  arrange(donor, metfiblabel)

frequency_sums <- frequency_table %>%
  group_by(donor) %>%
  summarise(total = sum(frequency), .groups = "drop")
if (any(abs(frequency_sums$total - 1) > 1e-12)) {
  stop("Donor-level frequencies do not sum to one.")
}

write.csv(
  transform(
    frequency_table,
    donor = as.character(donor),
    metfiblabel = as.character(metfiblabel),
    progression_status = as.character(progression_status)
  ),
  file.path(output_dir, "first_cohort_fibmet_frequencies_by_donor.csv"),
  row.names = FALSE
)

write.csv(
  transform(
    donor_annotation,
    donor = as.character(donor),
    progression_status = as.character(progression_status)
  ),
  file.path(output_dir, "first_cohort_donor_progression_annotation.csv"),
  row.names = FALSE
)

format_frequency <- function(x) {
  ifelse(
    x < 0.01,
    percent(x, accuracy = 0.1),
    percent(x, accuracy = 1)
  )
}

common_theme <- theme_classic(base_size = 15) +
  theme(
    axis.title = element_text(size = 16),
    axis.text.x = element_text(size = 13, face = "bold", angle = 35, hjust = 1),
    axis.text.y = element_text(size = 13),
    legend.position = "top",
    legend.box = "vertical",
    legend.justification = "center",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    plot.margin = margin(8, 16, 8, 8)
  )

status_layer <- geom_point(
  data = donor_annotation,
  aes(x = donor, y = 1.105, color = progression_status),
  inherit.aes = FALSE,
  shape = 15,
  size = 4.8
)

p_grouped <- ggplot(
  frequency_table,
  aes(x = donor, y = frequency, fill = metfiblabel)
) +
  geom_col(
    position = position_dodge(width = 0.78),
    width = 0.68,
    color = "grey25",
    linewidth = 0.35
  ) +
  geom_text(
    aes(label = format_frequency(frequency)),
    position = position_dodge(width = 0.78),
    vjust = -0.35,
    size = 4
  ) +
  status_layer +
  scale_fill_manual(
    name = "FibMET state",
    values = state_colors,
    breaks = state_levels
  ) +
  scale_color_manual(
    name = "Clinical course",
    values = status_colors,
    breaks = status_levels,
    labels = c("Progressive", "Stable")
  ) +
  scale_y_continuous(
    labels = label_percent(accuracy = 1),
    limits = c(0, 1.16),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0, 0))
  ) +
  guides(
    fill = guide_legend(order = 1, nrow = 1),
    color = guide_legend(order = 2, nrow = 1, override.aes = list(shape = 15, size = 5))
  ) +
  labs(
    title = "FibMET-state frequencies by donor",
    subtitle = "First cohort; colored squares indicate clinical course",
    x = NULL,
    y = "Frequency among stromal fibroblasts"
  ) +
  common_theme

stack_labels <- frequency_table %>%
  mutate(label = if_else(frequency >= 0.05, format_frequency(frequency), ""))

p_stacked <- ggplot(
  stack_labels,
  aes(x = donor, y = frequency, fill = metfiblabel)
) +
  geom_col(
    position = position_stack(reverse = TRUE),
    width = 0.72,
    color = "white",
    linewidth = 0.5
  ) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5, reverse = TRUE),
    color = "white",
    fontface = "bold",
    size = 4
  ) +
  status_layer +
  scale_fill_manual(
    name = "FibMET state",
    values = state_colors,
    breaks = state_levels
  ) +
  scale_color_manual(
    name = "Clinical course",
    values = status_colors,
    breaks = status_levels,
    labels = c("Progressive", "Stable")
  ) +
  scale_y_continuous(
    labels = label_percent(accuracy = 1),
    limits = c(0, 1.16),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0, 0))
  ) +
  guides(
    fill = guide_legend(order = 1, nrow = 1),
    color = guide_legend(order = 2, nrow = 1, override.aes = list(shape = 15, size = 5))
  ) +
  labs(
    title = "FibMET-state composition by donor",
    subtitle = "First cohort; colored squares indicate clinical course",
    x = NULL,
    y = "Frequency among stromal fibroblasts"
  ) +
  common_theme

save_plot <- function(plot, stem, width = 8.7, height = 5.5) {
  pdf(
    file.path(output_dir, paste0(stem, ".pdf")),
    width = width,
    height = height,
    useDingbats = FALSE
  )
  print(plot)
  dev.off()

  svg(
    file.path(output_dir, paste0(stem, ".svg")),
    width = width,
    height = height,
    bg = "white"
  )
  print(plot)
  dev.off()

  ggsave(
    file.path(output_dir, paste0(stem, ".png")),
    plot,
    width = width,
    height = height,
    dpi = 320,
    bg = "white"
  )
}

save_plot(p_grouped, "first_cohort_fibmet_frequency_grouped")
save_plot(p_stacked, "first_cohort_fibmet_frequency_stacked")

print(frequency_table)
