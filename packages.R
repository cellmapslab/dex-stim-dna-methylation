LoadPackages <- function(pkg.list){
    new.pkg <- pkg.list[!(pkg.list %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg.list, library, character.only = TRUE)
}

LoadPackagesBio <- function(pkg.list){
    new.pkg <- pkg.list[!(pkg.list %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
        BiocManager::install(new.pkg, dependencies = TRUE)
    sapply(pkg.list, library, character.only = TRUE)
}

pkg.list <- c("BiocManager", "future", "RColorBrewer", "RPMM", "ggplot2", "matrixStats", "reshape", "Hmisc", "factoextra")

biocmanager.pkg.list <- c("minfi", "missMethyl", "wateRmelon", "sva")

LoadPackages(pkg.list)
LoadPackagesBio(biocmanager.pkg.list)
