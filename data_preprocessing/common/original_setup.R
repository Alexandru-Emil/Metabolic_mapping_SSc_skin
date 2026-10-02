# Original IMC package/helper setup with explicit dependency checking.
source("data_preprocessing/common/config.R")
source("data_preprocessing/common/dependencies.R")
source("data_preprocessing/common/original_utils.R")
all_packages<-unique(c(required_packages_cran,required_packages_github,required_packages_bioconductor))
require_packages(all_packages)
invisible(lapply(all_packages,library,character.only=TRUE))
