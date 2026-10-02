# Dependency declarations from the original IMC config.R.
required_packages_cran <- c("renv", "here", "tidyverse", "knitr", "rmarkdown", "readr", "remotes", "gridExtra", 
    "grid", "pals", "cowplot", "RColorBrewer", "viridis", "rlang", "data.table", "chameleon", "circlize", 
    "fmsb", "ggstatsplot", "ggpubr", "rstatix", "ggrepel", "stringr", "spdep", "ggraph", "tidygraph", 
    "igraph", "Rtsne", "caret", "future", "future.apply", "progressr", "parallelly", "ggsci")
required_packages_github <- c("Seurat", "SeuratData", "SeuratWrappers", "Banksy", "ggradar", "patchwork", 
    "Rphenograph", "presto")
required_packages_bioconductor <- unique(c("SpatialExperiment", "SingleCellExperiment", "dittoSeq", "SpatialFeatureExperiment", 
    "scran", "SFEData", "Voyager", "scuttle", "scater", "ComplexHeatmap", "imcRtools", "CATALYST", "bluster", 
    "BiocParallel", "BiocSingular", "SingleR", "GenomeInfoDbData"))
