###UMAP and heatmaps##
library(RColorBrewer)
library(pals)

z<-list("glycolysis"=c(2,30,32,7,15,26),"TCA_cycle"=c(29,24,31,17),"Amino_Acid"=c(8),"Fatty_acid_synthesis"=c(36),
        "Fatty_acid_oxidation"=34,"Pentose_Pathway"=12,"anabolic_signaling"=c(9),
        "catabolic_signaling"=c(21),"Oxidative_stress"=c(37),
        "hypoxia"=38,"Mitochondrial_biogenesis"=c(22,13,10))


markers_ordered <- names(spe)[do.call(c, z)]
markers_ordered <- rownames(spe)[unlist(z, use.names = FALSE)]
markers_ordered

rowData(spe)$use_channel_met_markers <- grepl("Glut1|PKM2|CD98|p_mTOR|NRF2_p|G6PD|TFAM|LDHA|SDHA|pAMPK|PGC-1a|OGDH|GAPDH|CS|HK1|ATP5A|PFKL1|CPT1a|ACAC|NOX4|HIF1a", rownames(spe)) ## for old IMC
rowData(spe)$use_channel_met_markers <- grepl("GLUT1|PKM2|CD98|p_mTOR|pNRF2|G6PD|mtTFA|LDHA|
                                              SDHA|pAMPK|PGC1a|OGDH|GAPDH|CS|HK1|ATP5A|PFKL1_PFKM|
                                              CPT1a|ACAC|NOX4|HIF1a", rownames(spe)) ## for new spe

met_markers <- c('GLUT1','PKM2','CD98','p_mTOR','pNRF2','G6PD','mtTFA',
                 'LDHA','SDHA','pAMPK','PGC1a','OGDH','GAPDH','CS','HK1',
                 'ATP5A','PFKL1_PFKM','CPT1a','ACAC','NOX4','HIF1a') ## For new IMC data

celltype_mean <- aggregateAcrossCells(as(spe[rowData(spe)$use_channel_met_markers,], "SingleCellExperiment"),  
                                              ids = spe$"Rphenograph_1000", 
                                              statistics = "mean", 
                                              use.assay.type = "z-rescaled")
#jpeg("Stromal_Rpheno_1000.jpg")

spe_test <- celltype_mean
spe_test$Rphenograph_1000 <- as.character(spe_test$Rphenograph_1000)
dittoHeatmap(spe_test,
             assay = "z-rescaled", 
             genes = markers_ordered,
             cluster_cols = T,
             cluster_rows = F,
             scale = "none",
             heatmap.colors = c(viridis(50), rep("#FDE725FF", 1)),
             annot.by = c("Rphenograph_1000", 'ncells'))
