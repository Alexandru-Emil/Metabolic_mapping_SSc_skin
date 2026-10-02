
## Requires:
##   - genes exist in rownames(sfe)  
##   - colData(sfe)$IMC_Rphenograph_500 exists
## ============================================

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

## ----------------------------
## USER INPUTS (edit these)
## ----------------------------
# Objects assumed already in memory:
#   sfe : SpatialFeatureExperiment

genes_use_subset <- c(
  "CTHRC1","POSTN","FN1","COL1A1","COL1A2","COL8A1","SFRP4","SFRP2","FAP","LRRC15",
  "CCL19","CXCL12","APOE","PI16","LGR5","COL18A1","WIF1"
)

library(Seurat)

seurat_obj <- as.Seurat(
  new_sfe,
  counts = "counts",
  data   = "logcounts"
)

Idents(seurat_obj) <- "IMC_Rphenograph_500"

DotPlot(
  object = subset(seurat_obj, idents = levels(Idents(seurat_obj))),
  features = genes_use_subset,
  col.min = 0.5,
  col.max = 2,
  dot.scale = 5
) +
  coord_flip() +
  labs(x = NULL, y = NULL)

