# Scientific workflow derived from IMC_metabolic_analysis_AEM/R/Scripts/output_dir_svg_plots_fx.R
# Run from the repository root with METABOLIC_INPUT_DIR set to the external project data.
source("data_preprocessing/common/config.R")
# Create output directory to save plots as svg files
output_dir_svg_plots_fx <- function(ts) {
  # Conditional check: use knitr::current_input() only when knitting
  if (!is.null(knitr::current_input())) {
    # Get the current RMarkdown filename without extension
    md_filename <- tools::file_path_sans_ext(basename(knitr::current_input()))
  } else {
    # Fallback filename for interactive sessions
    md_filename <- "interactive_session"
  }
  
  # Create a directory to save the SVG files
  output_dir <- file.path(paste0(md_filename, "_svg_files_", ts))
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Return the directory path
  return(output_dir)
}