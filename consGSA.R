# Description:
#     Performs consensus GSA using eight common (and relatively fast) methods that are
#     included in Piano.
#     Runs on maximum number of CPU nodes. Gives resList as output, which can be
#     used to make a plot using consGSAplot.
#
# Input:      Pval      named vector of gene-level adjusted P-values
#             FC        named vector of gene-level log2 fold-changes
#             gsc       gene-set, as loaded using Piano
#             nPerm     the number of permutations to use for gene sampling, default = 1000
#             gsSizeDn  cutoff-value for minimum number of genes in gene-sets, default = 5
#             gsSizeUp  cutoff-value for maximum number of genes in gene-sets, default = 500
#
# Output:     resList   results list
#
# 2018-03-13  Eduard Kerkhoven



consGSA <-
  function(Pval,
           FC,
           gsc,
           nPerm = 1000,
           gsSizeDn = 5,
           gsSizeUp = 500) {
    require(piano)
    require(parallel)
    require(snowfall)
    require(tidyr)
    
    # Find out number of processor cores, run GSA on all cores.
    cores <- as.numeric(detectCores())
    cat("Run GSA on", cores, "CPU cores.\n")

    cat("Running GSA 1/10 (mean)")
    gsaRes1 <- runGSA(
      Pval,
      FC,
      geneSetStat = "mean",
      gsc = gsc,
      nPerm = round(nPerm / cores) * cores,
      gsSizeLim = c(gsSizeDn, gsSizeUp),
      ncpus = cores
    )
    cat("Running GSA 2/10 (median)")
    gsaRes2 <- runGSA(
      Pval,
      FC,
      geneSetStat = "median",
      gsc = gsc,
      nPerm = round(nPerm / cores) * cores,
      gsSizeLim = c(gsSizeDn, gsSizeUp),
      ncpus = cores
    )
    cat("Running GSA 3/10 (sum)")
    gsaRes3 <- runGSA(
      Pval,
      FC,
      geneSetStat = "sum",
      gsc = gsc,
      nPerm = round(nPerm / cores) * cores,
      gsSizeLim = c(gsSizeDn, gsSizeUp),
      ncpus = cores
    )
    cat("Running GSA 4/10 (stouffer)")
    gsaRes4 <- runGSA(
      Pval,
      FC,
      geneSetStat = "stouffer",
      gsc = gsc,
      nPerm = round(nPerm / cores) * cores,
      gsSizeLim = c(gsSizeDn, gsSizeUp),
      ncpus = cores
    )
    cat("Running GSA 5/10 (tailStrength)")
    gsaRes5 <- runGSA(
      Pval,
      FC,
      geneSetStat = "tailStrength",
      gsc = gsc,
      nPerm = round(nPerm / cores) * cores,
      gsSizeLim = c(gsSizeDn, gsSizeUp),
      ncpus = cores
    )
    cat("Running GSA 6/10 (gsea)")
    gsaRes6 <- runGSA(
      FC,
      geneSetStat = "gsea",
      gsc = gsc,
      nPerm = round(nPerm / cores) * cores,
      gsSizeLim = c(gsSizeDn, gsSizeUp),
      ncpus = cores
    )
    cat("Running GSA 7/10 (fisher)")
    gsaRes7 <- runGSA(
      Pval,
      FC,
      geneSetStat = "fisher",
      gsc = gsc,
      nPerm = round(nPerm / cores) * cores,
      gsSizeLim = c(gsSizeDn, gsSizeUp),
      ncpus = cores
    )
    cat("Running GSA 8/10 (maxmean)")
    gsaRes8 <- runGSA(
      FC,
      geneSetStat = "maxmean",
      gsc = gsc,
      nPerm = round(nPerm / cores) * cores,
      gsSizeLim = c(gsSizeDn, gsSizeUp),
      ncpus = cores
    )
    cat("Running GSA 9/10 (fgsea)")
    gsaRes9 <- runGSA(
      FC,
      geneSetStat = "fgsea",
      gsc = gsc,
      nPerm = round(nPerm / cores) * cores,
      gsSizeLim = c(gsSizeDn, gsSizeUp),
      ncpus = cores
    )
        cat("Running GSA 10/10 (page)")
    gsaRes10 <- runGSA(
      FC,
      geneSetStat = "page",
      gsc = gsc,
      nPerm = round(nPerm / cores) * cores,
      gsSizeLim = c(gsSizeDn, gsSizeUp),
      ncpus = cores
    )
    # No maxmean and fisher: don't support direction. No wilcoxon, too slow.
    cat("Reorganizing data and prepare for plotting")
    # Combine results in list
    resList <- list(gsaRes1, gsaRes2, gsaRes3, gsaRes4, gsaRes5, gsaRes6, gsaRes7,
                    gsaRes8, gsaRes9, gsaRes10)
    resList <-
      setNames(resList,
               c("mean", "median", "sum", "stouffer",
                 "tailStrength", "gsea", "fisher", "maxmean", "fgsea", "page"))
    
return(resList)
}