# This script installs packages that I routinely use.
# Run this after a clean installation of R.
# 2018-03-13  Eduard Kerkhoven
source("https://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")
biocLite(c("piano","limma","edgeR","ggplot2","snowfall","plyr","dplyr","tidyr","gdata","affy","plier","DESeq"))
install.packages(c("labeling","devtools","VennDiagram"))