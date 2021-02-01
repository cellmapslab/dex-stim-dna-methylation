LoadPackages <- function(pkg.list){
    new.pkg <- pkg.list[!(pkg.list %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
        install.packages(new.pkg, dependencies = TRUE)
    suppressMessages(sapply(pkg.list, library, character.only = TRUE))
}

LoadPackagesBio <- function(pkg.list){
    new.pkg <- pkg.list[!(pkg.list %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
        BiocManager::install(new.pkg, dependencies = TRUE)
    suppressMessages(sapply(pkg.list, library, character.only = TRUE))
}

pkg.list <- c("BiocManager", "future", "RColorBrewer", "RPMM", "ggplot2", "gplots", "matrixStats", 
              "reshape", "Hmisc", "factoextra", "ggpubr", "gridExtra", "grid", "tidyverse")

biocmanager.pkg.list <- c("minfi", "missMethyl", "wateRmelon", "sva", "mixOmics", "ENmix", "FlowSorted.Blood.EPIC", "FlowSorted.Blood.450k", "SummarizedExperiment")

LoadPackages(pkg.list)
LoadPackagesBio(biocmanager.pkg.list)
