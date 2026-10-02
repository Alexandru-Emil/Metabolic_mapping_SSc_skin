# Paths are configured by the caller. No packages are installed by this file.
project_dir <- Sys.getenv("METABOLIC_INPUT_DIR", unset = "inputs")
analysis_output_dir <- Sys.getenv("METABOLIC_OUTPUT_DIR", unset = "analysis_outputs")
project_path <- function(relative = "") {
  if (grepl("\\.(R|Rmd)$", relative, ignore.case = TRUE)) {
    name <- gsub(" ", "_", basename(relative), fixed = TRUE)
    candidates <- list.files("data_preprocessing", recursive = TRUE, full.names = TRUE)
    candidate_names<-sub("^[a-f0-9]{8}_","",basename(candidates))
    hit <- candidates[candidate_names == name]
    if (length(hit) == 1L) return(hit)
  }
  if(grepl("R/Scripts/?$",gsub("\\\\","/",relative)))return("data_preprocessing/common")
  file.path(project_dir, relative)
}
dir.create(analysis_output_dir, recursive = TRUE, showWarnings = FALSE)
require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install the documented dependencies before running: ", paste(missing, collapse = ", "))
}
read_analysis_object <- function(path, object = "spe") {
  if (!file.exists(path)) stop("Input not found: ", path)
  if (grepl("\\.RData$", path, ignore.case = TRUE)) {
    env <- new.env(parent = baseenv()); load(path, envir = env)
    if (exists(object, envir = env, inherits = FALSE)) return(env[[object]])
    if (exists("spe_list", envir = env, inherits = FALSE)) return(env$spe_list$spe)
    stop("Object ", object, " not found in workspace")
  }
  readRDS(path)
}
validation_donors <- paste0("Validation", 1:6)
first_cohort_donors <- c(paste0("Healthy", 1:7), paste0("SSc", 1:9))
metabolic_markers_reference <- c("Glut1", "PKM2", "CD98", "p-mTOR", "NRF2-p", "G6PD", "TFAM", "LDHA", "SDHA", "pAMPK", "PGC-1a", "OGDH", "GAPDH", "CS", "HK1", "ATP5A", "PFKL1", "CPT1a", "ACAC", "NOX4", "HIF1a")
glycolysis_markers_reference <- c("Glut1", "HK1", "PFKL1", "PKM2", "LDHA", "GAPDH")
tca_oxphos_markers_reference <- c("CS", "OGDH", "ATP5A", "SDHA")
