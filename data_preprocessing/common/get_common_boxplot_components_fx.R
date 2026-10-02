# Scientific workflow derived from IMC_metabolic_analysis_AEM/R/Scripts/get_common_boxplot_components_fx.R
# Run from the repository root with METABOLIC_INPUT_DIR set to the external project data.
source("data_preprocessing/common/config.R")
# Define a function to return common plot components to use with boxplots
get_common_boxplot_components <- function() {
  list(
    geom_point(position = position_dodge2(width = 0.5, padding = 1.5)),
    theme_minimal(),
    theme(
      axis.text.x = element_text(size = 14, face = "bold", angle = 45, hjust = 1),
      axis.text.y = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      legend.title = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 16, face = "bold"),
      strip.text = element_text(size = 14, face = "bold")
    )
  )
}