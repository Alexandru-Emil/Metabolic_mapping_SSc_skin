# Scientific workflow derived from Cytomapper R scripts Veda, Nan/z score rescaling function.R
# Run from the repository root with METABOLIC_INPUT_DIR set to the external project data.
source("data_preprocessing/common/config.R")
do.rescale.z <- function(dat,
                       use.cols,
                       append.name = '_rescaled_z'
){ 
  
  ### Packages
  require('data.table')
  
  ### Test data
  # dat <- Spectre::demo.asinh
  # use.cols <- names(dat)[c(11:19)]
  # new.min = -1
  # new.max = 1
  # append.name = '_rescaled'
  
  ### Create normalisation function
  #norm.fun <- function(x) {(x - min(x, na.rm=TRUE))/(max(x,na.rm=TRUE) -min(x, na.rm=TRUE))}
  norm.fun <- function(x) {(x - mean(x))/(sd(x))}
  #norm.fun <- function(x){scales::rescale(value, to=c(new.min,new.max))}
  
  ### Establish dataset
  value <- dat[,use.cols,with = FALSE]
  
  ### Normalise between new values
  res <- as.data.table(lapply(value, norm.fun)) # by default, removes the names of each row
  names(res) <- paste0(names(res), append.name)
  dat <- cbind(dat, res)
  
  ### Return
  return(dat)
}