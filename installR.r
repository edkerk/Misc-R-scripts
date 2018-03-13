# This script installs packages that I routinely use.
# Run this after a clean installation of R.
# 2018-03-13  Eduard Kerkhoven
.libPaths("C:/Work/Rlib")
source("https://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")
biocLite(c("piano","limma","edgeR","snowfall","tidyverse","gdata","affy","plier","DESeq"))
install.packages(c("labeling","devtools","VennDiagram","cluster","lattice","survival"))
