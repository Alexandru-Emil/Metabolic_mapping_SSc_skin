# Scientific workflow derived from Cytomapper R scripts Veda, Nan/Rescaling in sce or spe.R
# Run from the repository root with METABOLIC_INPUT_DIR set to the external project data.
source("data_preprocessing/common/config.R")
# Min-max rescaling
assay(spe, "rescaled") <- t(apply(assay(spe, "asinh"), 1, function(x){
  (x - min(x))/(max(x)-min(x))}))

## Experimental !! Z-score rescaling
assay(spe, "rescaled.z") <- t(apply(assay(spe, "asinh"), 1, function(x){
  (x - mean(x))/(sd(x))}))
