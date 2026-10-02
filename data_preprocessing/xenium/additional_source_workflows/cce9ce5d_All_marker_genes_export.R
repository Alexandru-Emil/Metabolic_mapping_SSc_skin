source("data_preprocessing/common/config.R")
data_dir <- '/Volumes/T7/Dataset/Level 2/Repeat10_Stromal/Re0.5 lamda0.3' 

input_spe_filename <- "sfe_merged derma_20250510180909_20250602133429.rds"  

spe <- readRDS(here(data_dir, input_spe_filename))

marker_table <- as.data.frame(metadata(spe)[["BANKSY_params_level_2_repeat_10"]][["markers_seurat"]])

write_csv(marker_table, project_path("external/lv2_repeat10_stromal_Re0.5 lam0.3.csv"))



