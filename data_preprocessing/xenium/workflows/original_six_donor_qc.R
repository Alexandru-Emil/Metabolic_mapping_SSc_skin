# Executed R code recovered from the original sfe_QC.html analysis.
source("data_preprocessing/common/config.R")
ptm <- proc.time()
tstamp <- format(Sys.time(), '%y%m%d%H%M%S')

library(SpatialFeatureExperiment)
library(tidyverse)

library(dplyr)
library(scales)
library(viridis)
library(SingleCellExperiment)

sfe <- readr::read_rds(project_path("Revision/Xenium/sfe object/sfe object250723154452.rds"))

counts_mat <- assays(sfe)[["counts"]]
nFeature <- colSums(counts_mat > 0)
colData(sfe)$nFeature <- nFeature


colData(sfe)$new_prop_neg_count <- 1 - (colData(sfe)$transcript_counts / 
                                          (colData(sfe)$control_probe_counts +
                                             colData(sfe)$genomic_control_counts +
                                             colData(sfe)$control_codeword_counts +
                                             colData(sfe)$unassigned_codeword_counts +
                                             colData(sfe)$transcript_counts))

# seurat$SampleId <- as.factor(seurat$SampleId)
library(ggplot2)
library(ggforce)


df <- as.data.frame(colData(sfe))



ggplot(df, aes(x = SampleId, y = transcript_counts)) +
  geom_violin(fill = "lightblue", alpha = 1) +
  geom_jitter(width = 0.4, size = 0.1, alpha = 0.2) +  # 模拟 pt.size 和 alpha
  geom_hline(yintercept = 20, linetype = "dashed", color = "red", linewidth = 1) +
  scale_y_continuous(trans = "pseudo_log") +
  theme_bw() +
  labs(title = "Transcript counts per cell", y = "Transcript Counts", x = "Sample ID")

# ggplot(seurat@meta.data, aes(x = SampleId, y = nCount_Xenium, fill = SampleId)) + 
#   geom_violin() + 
#   geom_sina(alpha = .3, size = .001, color = 'black') + 
#   geom_hline(yintercept = 20, color = 'red', linetype)

ggplot(df, aes(x = transcript_counts)) +
  geom_histogram(binwidth = 1, fill = "burlywood3") +
  scale_x_continuous(limits = c(0, 1000)) +
  theme_bw() +
  geom_vline(xintercept = 20, color = "red", linetype = "dashed", size = 1) 

ggplot(df, aes(x = SampleId, y = nFeature)) +
  geom_violin(fill = "lightblue", alpha = 1) +
  geom_jitter(width = 0.4, size = 0.1, alpha = 0.2) +  # 模拟 pt.size 和 alpha
  geom_hline(yintercept = 10, linetype = "dashed", color = "red", linewidth = 1) +
  scale_y_continuous(trans = "pseudo_log") +
  theme_bw() +
  labs(title = "nFeature", y = "Feature Counts", x = "Sample ID")

ggplot(df, aes(x = nFeature)) +
  geom_histogram(binwidth = 1, fill = "burlywood3") +
  scale_x_continuous(limits = c(0, 500)) +
  theme_bw() +
  geom_vline(xintercept = 10, color = "red", linetype = "dashed", size = 1) 

ggplot(df, aes(x = SampleId, y = cell_area)) +
  geom_violin(fill = "lightblue", alpha = 1) +
  geom_jitter(width = 0.4, size = 0.1, alpha = 0.2) +  
  geom_hline(yintercept = 300, linetype = "dashed", color = "red", linewidth = 1) +
  scale_y_continuous(trans = "pseudo_log") +
  theme_bw() +
  labs(title = "cell area", y = "cell area", x = "Sample ID")

ggplot(df, aes(x = SampleId, y = cell_area)) +
  geom_violin(fill = "lightblue", alpha = 1) +
  geom_jitter(width = 0.4, size = 0.1, alpha = 0.2) +  
  geom_hline(yintercept = 15, linetype = "dashed", color = "red", linewidth = 1) +
  scale_y_continuous(trans = "pseudo_log") +
  theme_bw() +
  labs(title = "cell area", y = "cell area", x = "Sample ID")

ggplot(df, aes(x = cell_area)) +
  geom_histogram(binwidth = 1, fill = "burlywood3") +
  scale_x_continuous(limits = c(0, 500)) +
  theme_bw() +
  geom_vline(xintercept = c(15, 300), color = "red", linetype = "dashed", size = 1) 

area_thres <- list(
  min_thres = c(Validation5 = 15,
                Validation6 = 15,
                'Validation1' = 19,
                'Validation2' = 13,
                Validation3 = 13,
                Validation4 = 14
),
  max_thres = c(Validation5 = 270,
                Validation6 = 250,
                'Validation1' = 220,
                'Validation2' = 300,
                Validation3 = 280,
                Validation4 = 300
  ))

sample_ids <- sfe$SampleId
area_vals <- sfe$cell_area

qc_pass <- mapply(function(id, area) {
  min <- area_thres$min_thres[[id]]
  max <- area_thres$max_thres[[id]]
  area >= min & area <= max
}, sample_ids, area_vals)

colData(sfe)$qc_pass_area <- qc_pass

ggplot(df, aes(x = SampleId, y = new_prop_neg_count)) +
  geom_violin(fill = "lightblue", alpha = 0.3, scale = "width") +
  geom_jitter(width = 0.3, size = 0.2, alpha = 0.5) +
  geom_hline(yintercept = 0.001, linetype = "dashed", color = "red", size = 1) +
   scale_y_continuous(trans = pseudo_log_trans(base = 10)) +
  coord_cartesian(ylim = c(1e-07, 5e-3)) +
  theme_bw() +
  labs(title = "prop_neg_count",
       y = "prop_neg_count (log10)", x = "Sample ID")

ggplot(df, aes(x = new_prop_neg_count)) +
  geom_histogram(binwidth = .001, fill = "burlywood3") +
  scale_x_continuous(limits = c(-.005, .02)) +
  theme_bw() +
  geom_vline(xintercept = c(.001), color = "red", linetype = "dashed", size = 1) 

sfe@colData$qc_pass <- case_when(sfe@colData@listData[["transcript_counts"]] > 20 &
                                        sfe@colData@listData[["nFeature"]] > 10 &
                                        sfe@colData@listData[["qc_pass_area"]] &
                                        sfe@colData@listData[["new_prop_neg_count"]] < .001 ~ TRUE, 
                                      .default = FALSE)

sfe$qc_pass <- factor(sfe$qc_pass, levels = c(TRUE, FALSE))  

samples <- unique(colData(sfe)$SampleId)
coords <- as.data.frame(spatialCoords(sfe))
coldata <- as.data.frame(colData(sfe))
coldata$spatialX <- coords$x_centroid
coldata$spatialY <- coords$y_centroid

for (i in samples) {
  cat("\n\n### Sample:", i, "\n")
  
  subdata <- coldata[coldata$SampleId == i, ]
  
  p1 <- ggplot(subdata, aes(x = spatialX, y = spatialY, color = qc_pass)) +
    geom_point(size = 0.3, alpha = 0.6) +  
    ggtitle(paste("QC Pass status for", i)) +
    theme_minimal()
  print(p1)  
  
p2 <- ggplot(subdata, aes(x = spatialX, y = spatialY, color = transcript_counts)) +
  geom_point(size = 0.3, alpha = 0.6) +
  scale_color_viridis_c(limits = c(min(subdata$transcript_counts), quantile(subdata$transcript_counts, 0.8))) +
  ggtitle(paste("transcript_counts in", i)) +
  theme_minimal()
 print(p2)
}

sfe <- sfe[, sfe@colData$qc_pass == TRUE]

saveRDS(sfe, paste0(project_path("Revision/Xenium/sfe object/sfe_afterqc"), tstamp, '.rds'))

etm <- proc.time() - ptm
cat('This analysis took :', etm[3], 's \n')

sessioninfo::session_info()