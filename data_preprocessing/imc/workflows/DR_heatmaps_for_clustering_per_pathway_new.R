# Scientific workflow derived from Older data and analysis/DR_heatmaps_for_clustering_per_pathway_new.R
# Run from the repository root with METABOLIC_INPUT_DIR set to the external project data.
source("data_preprocessing/common/config.R")
library(SpatialExperiment)
library(cytomapper)
library(HDF5Array)
library(dittoSeq)
library(uwot)
library(scater)
library(CATALYST)
library(Rphenograph)
library(igraph)
library(dplyr)

load(project_path("final_spe.RData"))

colvector <- colnames(final_spe)
spe_merged <- final_spe

rowData(spe_merged)$use_channel_met_markers <- grepl("Glut1|PKM2|CD98|p_mTOR|NRF2_p|G6PD|TFAM|LDHA|SDHA|pAMPK|PGC-1a|OGDH|GAPDH|CS|HK1|ATP5A|PFKL1|CPT1a|ACAC|NOX4|HIF1a", rownames(spe_merged))

z<-list("glycolysis"=c(2,30,32,7,15,26),"TCA_cycle"=c(29,24,31,17),"Amino_Acid"=c(8),"Fatty_acid_synthesis"=c(36),
        "Fatty_acid_oxidation"=34,"Pentose_Pathway"=12,"anabolic_signaling"=c(9),
        "catabolic_signaling"=c(21),"Oxidative_stress"=c(37),
        "hypoxia"=38,"Mitochondrial_biogenesis"=c(22,13,10))

reducedDim(spe_merged,"testUMAP") <- matrix(rep(NA,2*ncol(spe_merged)), nrow = ncol(spe_merged), ncol = 2)

for(i in unique(spe_merged$cell_labels)){
  w<-which(spe_merged$cell_labels==i)
  spe_sub<-spe_merged[,w]
  
  mat_x <- assay(spe_sub,"z_c_rescaled") %>% as.data.frame()
  
  for(j in 1:length(z)){
    if(length(z[[j]])==1){
      mat_x[names(z[j]),]<-assay(spe_sub,"z_c_rescaled")[z[[j]],]
    }
    else{
      mat_x[names(z[j]),]<-colMeans(assays(spe_sub)[[5]][z[[j]],])
    }
  }
  mat_x <- as.matrix(mat_x)
  UMAP_x <- calculateUMAP(mat_x, subset_row = names(z))
  reducedDim(spe_sub, "testUMAP") <- UMAP_x

  reducedDim(spe_merged,"testUMAP")[w,]<-reducedDim(spe_sub, "testUMAP")
  cat(i)
}

spe_stromal <- spe_merged[,spe_merged$cell_labels == "Stromal"]
spe_CD31 <- spe_merged[,spe_merged$cell_labels == "CD31+"]
spe_CD45 <- spe_merged[,spe_merged$cell_labels == "CD45+"]
spe_Ecad <- spe_merged[,spe_merged$cell_labels == "E-cad+"]
spe_unlabeled <- spe_merged[,spe_merged$cell_labels == "unlabeled"]


setwd(project_path("test plots"))

spe_stromal$Rphenograph_cluster_1000<-as.character(spe_stromal$Rphenograph_cluster_1000)
spe_CD31$Rphenograph_cluster_1000<-as.character(spe_CD31$Rphenograph_cluster_1000)
spe_CD45$Rphenograph_cluster_1000<-as.character(spe_CD45$Rphenograph_cluster_1000)
spe_Ecad$Rphenograph_cluster_1000<-as.character(spe_Ecad$Rphenograph_cluster_1000)

jpeg("Stromal_UMAP_Met.jpg")
dittoDimPlot(spe_stromal, var = "Rphenograph_cluster_1000", 
             reduction.use = "testUMAP", size = 0.2,
             do.label = TRUE) +
  ggtitle("Clusters expression on UMAP Met stromal")
dev.off()

jpeg("CD31_UMAP_Met.jpg")
dittoDimPlot(spe_CD31, var = "Rphenograph_cluster_1000", 
             reduction.use = "testUMAP", size = 0.2,
             do.label = TRUE) +
  ggtitle("Clusters expression on UMAP Met stromal")
dev.off()

jpeg("CD45_UMAP_Met.jpg")
dittoDimPlot(spe_CD45, var = "Rphenograph_cluster_1000", 
             reduction.use = "testUMAP", size = 0.2,
             do.label = TRUE) +
  ggtitle("Clusters expression on UMAP Met stromal")
dev.off()

jpeg("Ecad_UMAP_Met.jpg")
dittoDimPlot(spe_Ecad, var = "Rphenograph_cluster_1000", 
             reduction.use = "testUMAP", size = 0.2,
             do.label = TRUE) +
  ggtitle("Clusters expression on UMAP Met stromal")
dev.off()

spe_merged_v1 <- spe_merged

library(RColorBrewer)
library(pals)

pal.bands(alphabet, alphabet2, cols25, glasbey, kelly, polychrome, 
          stepped, tol, watlington,
          show.names=FALSE)

patient_id <- setNames(cols25(n= length(unique(spe_merged_v1$patient_id))), 
                       unique(spe_merged_v1$patient_id))

sample_id <- setNames(c(brewer.pal(6, "YlOrRd")[3:6],brewer.pal(6, "PuBu")[3:6],brewer.pal(6, "YlGn")[3:6],brewer.pal(6, "BuPu")[3:6]),unique(spe_merged_v1$sample_id))

indication <- setNames(c("#8DD3C7", "#FFFFB3"), unique(spe_merged_v1$indication))

color_vectors <- list()

color_vectors$patient_id <- patient_id
color_vectors$sample_id <- sample_id
color_vectors$indication <- indication

metadata(spe_stromal)$color_vectors <- color_vectors
metadata(spe_CD31)$color_vectors <- color_vectors
metadata(spe_CD45)$color_vectors <- color_vectors
metadata(spe_Ecad)$color_vectors <- color_vectors

celltype <- setNames(c("#3F1B03", "#F4AD31", "#894F36", "#1C750C"),
                     c("Stromal", "CD31+", "CD45+", "E-cad+"))

metadata(spe_stromal)$color_vectors$celltype <- celltype
metadata(spe_CD31)$color_vectors$celltype <- celltype
metadata(spe_CD45)$color_vectors$celltype <- celltype
metadata(spe_Ecad)$color_vectors$celltype <- celltype

markers_ordered <- names(spe_merged)[do.call(c, z)]


## Stromal ##
celltype_mean_stromal <- aggregateAcrossCells(as(spe_stromal[rowData(spe_stromal)$use_channel_met_markers,], "SingleCellExperiment"),  
                                      ids = spe_stromal$"Rphenograph_cluster_1000", 
                                      statistics = "mean", 
                                      use.assay.type = "z_c_rescaled")
jpeg("Stromal_Rpheno_1000.jpg")
dittoHeatmap(celltype_mean_stromal,
             assay = "z_c_rescaled", 
             genes = markers_ordered,
             cluster_cols = TRUE,
             cluster_rows = F,
             scale = "none",
             heatmap.colors = viridis(100),
             annot.by = c("Rphenograph_cluster_1000", "ncells"),
             annot.colors = c(dittoColors(1)[1:length(unique(spe_stromal$"Rphenograph_cluster_1000"))], dittoColors(1)[5]))
dev.off()

## CD31 ##
celltype_mean_CD31 <- aggregateAcrossCells(as(spe_CD31[rowData(spe_CD31)$use_channel_met_markers,], "SingleCellExperiment"),  
                                      ids = spe_CD31$"Rphenograph_cluster_1000", 
                                      statistics = "mean", 
                                      use.assay.type = "z_c_rescaled")
jpeg("CD31_Rpheno_1000.jpg")
dittoHeatmap(celltype_mean_CD31,
             assay = "z_c_rescaled",
             genes = markers_ordered,
             cluster_cols = TRUE, 
             cluster_rows = F,
             scale = "none",
             heatmap.colors = viridis(100),
             annot.by = c("Rphenograph_cluster_1000", "ncells"),
             annot.colors = c(dittoColors(1)[1:length(unique(spe_CD31$"Rphenograph_cluster_1000"))], dittoColors(1)[5]))
dev.off()

## CD45 ##
celltype_mean_CD45 <- aggregateAcrossCells(as(spe_CD45[rowData(spe_CD45)$use_channel_met_markers,], "SingleCellExperiment"),  
                                      ids = spe_CD45$"Rphenograph_cluster_1000", 
                                      statistics = "mean", 
                                      use.assay.type = "z_c_rescaled")

jpeg("CD45_Rpheno_1000.jpg")
dittoHeatmap(celltype_mean_CD45,
             assay = "z_c_rescaled",
             genes = markers_ordered,
             cluster_cols = TRUE, 
             cluster_rows = F,
             scale = "none",
             heatmap.colors = viridis(100),
             annot.by = c("Rphenograph_cluster_1000", "ncells"),
             annot.colors = c(dittoColors(1)[1:length(unique(spe_CD45$"Rphenograph_cluster_1000"))], dittoColors(1)[5]))
dev.off()

## Ecad ##
celltype_mean_Ecad <- aggregateAcrossCells(as(spe_Ecad[rowData(spe_Ecad)$use_channel_met_markers,], "SingleCellExperiment"),  
                                      ids = spe_Ecad$"Rphenograph_cluster_1000", 
                                      statistics = "mean", 
                                      use.assay.type = "z_c_rescaled")
jpeg("Ecad_Rpheno_1000.jpg")
dittoHeatmap(celltype_mean_Ecad,
             assay = "z_c_rescaled",
             genes = markers_ordered,
             cluster_cols = TRUE, 
             cluster_rows = F,
             scale = "none",
             heatmap.colors = viridis(100),
             annot.by = c("Rphenograph_cluster_1000", "ncells"),
             annot.colors = c(dittoColors(1)[1:length(unique(spe_Ecad$"Rphenograph_cluster_1000"))], dittoColors(1)[5]))
dev.off()
