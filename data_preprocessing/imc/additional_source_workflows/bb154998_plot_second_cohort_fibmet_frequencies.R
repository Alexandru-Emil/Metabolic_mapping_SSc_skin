source("data_preprocessing/common/config.R")
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

spe_file <- paste0(
  paste0(project_path("Revision/"), "/"),
  "Nan annotation/Data/IMC spe objects/New IMC/spe2025-12-17 umap.RData"
)
output_dir <- file.path(getwd(), "second_cohort_fibmet_frequency_outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

load(spe_file)

sample_to_donor <- c(
  "Validation1" = "Validation1",
  "Validation2" = "Validation2",
  "ExcludedValidation" = "ExcludedValidation",
  "Validation3" = "Validation3",
  "Validation4" = "Validation4",
  "Validation5" = "Validation5",
  "Validation6" = "Validation6"
)

donor_levels <- c("Validation1", "Validation2", "ExcludedValidation", "Validation3", "Validation4", "Validation5", "Validation6")
state_levels <- c("Met_hi_Fib", "Other_Fib")
state_colors <- c("Met_hi_Fib" = "#E15759", "Other_Fib" = "#4E79A7")

metadata <- as.data.frame(colData(spe))
unmapped_samples <- setdiff(unique(as.character(metadata$sample_id)), names(sample_to_donor))
if (length(unmapped_samples) > 0) {
  stop("Missing donor labels for: ", paste(unmapped_samples, collapse = ", "))
}

fibmet_cells <- metadata %>%
  transmute(
    sample_id = as.character(sample_id),
    donor = unname(sample_to_donor[as.character(sample_id)]),
    gate = as.character(gates),
    metfiblabel = case_when(
      as.character(Rphenograph_500) %in% c("1", "2") ~ "Met_hi_Fib",
      as.character(Rphenograph_500) == "3" ~ "Other_Fib",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(metfiblabel))

if (!all(fibmet_cells$gate == "Stromal")) {
  stop("FibMET-labelled cells were found outside the Stromal gate.")
}

frequency_table <- fibmet_cells %>%
  count(donor, metfiblabel, name = "n_cells") %>%
  group_by(donor) %>%
  mutate(
    n_fibmet_cells = sum(n_cells),
    frequency = n_cells / n_fibmet_cells
  ) %>%
  ungroup() %>%
  mutate(
    donor = factor(donor, levels = donor_levels),
    metfiblabel = factor(metfiblabel, levels = state_levels)
  ) %>%
  arrange(donor, metfiblabel)

if (any(abs(
  frequency_table %>%
    group_by(donor) %>%
    summarise(total = sum(frequency), .groups = "drop") %>%
    pull(total) - 1
) > 1e-12)) {
  stop("Donor-level frequencies do not sum to one.")
}

write.csv(
  transform(
    frequency_table,
    donor = as.character(donor),
    metfiblabel = as.character(metfiblabel)
  ),
  file.path(output_dir, "second_cohort_fibmet_frequencies_by_donor.csv"),
  row.names = FALSE
)

common_theme <- theme_classic(base_size = 15) +
  theme(
    axis.title = element_text(size = 16),
    axis.text.x = element_text(size = 13, face = "bold", angle = 35, hjust = 1),
    axis.text.y = element_text(size = 13),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 13),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    plot.margin = margin(8, 16, 8, 8)
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
    aes(label = percent(frequency, accuracy = 1)),
    position = position_dodge(width = 0.78),
    vjust = -0.35,
    size = 4.2
  ) +
  scale_fill_manual(values = state_colors, breaks = state_levels) +
  scale_y_continuous(
    labels = label_percent(accuracy = 1),
    limits = c(0, 1.08),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "FibMET-state frequencies by donor",
    subtitle = "Second cohort",
    x = NULL,
    y = "Frequency among FibMET-labelled stromal cells"
  ) +
  common_theme

p_stacked <- ggplot(
  frequency_table,
  aes(x = donor, y = frequency, fill = metfiblabel)
) +
  geom_col(
    position = position_stack(reverse = TRUE),
    width = 0.72,
    color = "white",
    linewidth = 0.5
  ) +
  geom_text(
    aes(label = percent(frequency, accuracy = 1)),
    position = position_stack(vjust = 0.5, reverse = TRUE),
    color = "white",
    fontface = "bold",
    size = 4.2
  ) +
  scale_fill_manual(values = state_colors, breaks = state_levels) +
  scale_y_continuous(
    labels = label_percent(accuracy = 1),
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "FibMET-state composition by donor",
    subtitle = "Second cohort",
    x = NULL,
    y = "Frequency among FibMET-labelled stromal cells"
  ) +
  common_theme

save_plot <- function(plot, stem, width = 7.6, height = 4.9) {
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

save_plot(p_grouped, "second_cohort_fibmet_frequency_grouped")
save_plot(p_stacked, "second_cohort_fibmet_frequency_stacked")

print(frequency_table)
