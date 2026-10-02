# Scientific workflow derived from IMC_metabolic_analysis_AEM/R/Scripts/get_ts_rightformat_fx.R
# Run from the repository root with METABOLIC_INPUT_DIR set to the external project data.
source("data_preprocessing/common/config.R")
get_ts_rightformat <- function() {
  return(format(Sys.time(), "%Y%m%d%H%M%S"))
}