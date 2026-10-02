source("data_preprocessing/common/config.R")
library(ggplot2)
library(dplyr)

output_dir <- project_path("Revision/Nan annotation/Data/RMD files/scmetabolism_outputs")

fib_summary <- read.csv(file.path(output_dir, "scmetabolism_metfiblabel_energy_summary.csv"), check.names = FALSE)
ec_summary <- read.csv(file.path(output_dir, "scmetabolism_ec_metEClabel_energy_summary.csv"), check.names = FALSE)

pathway_levels <- c("Glycolysis / Gluconeogenesis", "Citrate cycle (TCA cycle)")
pathway_labels <- c(
  "Glycolysis / Gluconeogenesis" = "Glycolysis /\nGluconeogenesis",
  "Citrate cycle (TCA cycle)" = "Citrate cycle\n(TCA cycle)"
)

make_energy_dotplot_small <- function(
  summary_df,
  group_col,
  group_levels,
  title,
  scale_by_column = FALSE
) {
  plot_df <- summary_df %>%
    dplyr::filter(pathway %in% pathway_levels) %>%
    dplyr::mutate(
      pathway = factor(as.character(pathway), levels = pathway_levels),
      group_label = factor(as.character(.data[[group_col]]), levels = rev(group_levels))
    )

  if (isTRUE(scale_by_column)) {
    plot_df <- plot_df %>%
      dplyr::group_by(pathway) %>%
      dplyr::mutate(
        dot_min = min(mean_score, na.rm = TRUE),
        dot_max = max(mean_score, na.rm = TRUE),
        dot_score = dplyr::if_else(
          dot_max > dot_min,
          (mean_score - dot_min) / (dot_max - dot_min),
          0.5
        )
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(-dot_min, -dot_max)
    legend_title <- "Column-scaled mean"
  } else {
    plot_df <- plot_df %>% dplyr::mutate(dot_score = mean_score)
    legend_title <- "Mean AUCell score"
  }

  ggplot(plot_df, aes(x = pathway, y = group_label)) +
    geom_point(
      aes(fill = dot_score, size = dot_score),
      shape = 21,
      color = "grey35",
      stroke = 0.35
    ) +
    scale_x_discrete(labels = pathway_labels) +
    scale_fill_gradient(low = "white", high = "#D7301F", name = legend_title) +
    scale_size(range = c(5.8, 11.5), guide = "none") +
    labs(title = title, x = NULL, y = NULL) +
    theme_bw(base_size = 13) +
    theme(
      axis.text.x = element_text(size = 12, lineheight = 0.95),
      axis.text.y = element_text(size = 13, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10.5),
      panel.grid.major = element_line(color = "grey88"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 15, hjust = 0.5, margin = margin(b = 6)),
      plot.margin = margin(5.5, 6, 5.5, 5.5)
    )
}

save_vector_plot <- function(plot, basename, width = 4.9, height = 2.65) {
  pdf_file <- file.path(output_dir, paste0(basename, ".pdf"))
  svg_file <- file.path(output_dir, paste0(basename, ".svg"))

  grDevices::pdf(pdf_file, width = width, height = height, useDingbats = FALSE)
  print(plot)
  grDevices::dev.off()

  grDevices::svg(svg_file, width = width, height = height)
  print(plot)
  grDevices::dev.off()

  invisible(c(pdf_file, svg_file))
}

plots <- list(
  scmetabolism_metfiblabel_energy_dotplot_raw_small = make_energy_dotplot_small(
    fib_summary,
    group_col = "metfiblabel",
    group_levels = c("Other_Fib", "Met_hi_Fib"),
    title = "Fibroblast mean Glycolysis and TCA scores",
    scale_by_column = FALSE
  ),
  scmetabolism_metfiblabel_energy_dotplot_scaled_by_column_small = make_energy_dotplot_small(
    fib_summary,
    group_col = "metfiblabel",
    group_levels = c("Other_Fib", "Met_hi_Fib"),
    title = "Fibroblast mean Glycolysis and TCA scores - scaled by column",
    scale_by_column = TRUE
  ),
  scmetabolism_ec_metEClabel_energy_dotplot_raw_small = make_energy_dotplot_small(
    ec_summary,
    group_col = "metEClabel",
    group_levels = c("Met_int_EC", "Met_hi_EC"),
    title = "EC mean Glycolysis and TCA scores",
    scale_by_column = FALSE
  ),
  scmetabolism_ec_metEClabel_energy_dotplot_scaled_by_column_small = make_energy_dotplot_small(
    ec_summary,
    group_col = "metEClabel",
    group_levels = c("Met_int_EC", "Met_hi_EC"),
    title = "EC mean Glycolysis and TCA scores - scaled by column",
    scale_by_column = TRUE
  )
)

generated <- unlist(Map(save_vector_plot, plots, names(plots)))
writeLines(generated)
