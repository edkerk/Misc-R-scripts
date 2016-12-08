# Description:
#     Performs consensus GSA using five common (and relatively fast) methods that are
#     included in Piano: "mean", "median", "sum", "stouffer" and "tailStrength".
#     Runs on maximum number of CPU nodes. Gives resList as output, which can be
#     used to make a plot using consGSAplot.
#
# Input:      Pval      named vector of gene-level adjusted P-values
#             FC        named vector of gene-level log2 fold-changes
#             GS        gene-set, as loaded usign Piano
# Output:     resList   results list
#
# 2016-11-28 EJK: Separate funtion to run consensus GSA.
# 2016-05-16 Eduard Kerkhoven (eduardk@chalmers.se)



consGSA <-
  function(Pval,
           FC,
           GS) {
    require(piano)
    require(parallel)
    require(snowfall)
    require(tidyr)
    #  if (!exists("rankScore")) # Attempt to set default settings, not sure how to do this...
    
    # Find out number of processor cores, run GSA on n-1 cores.
    cores <- as.numeric(detectCores())
    cat("Running GSA 1/5")
    gsaRes1 <- runGSA(
      Pval,
      FC,
      geneSetStat = "mean",
      gsc = GS,
      nPerm = round(1000 / cores) * cores,
      gsSizeLim = c(5, 200),
      ncpus = cores
    )
    cat("Running GSA 2/5")
    gsaRes2 <- runGSA(
      Pval,
      FC,
      geneSetStat = "median",
      gsc = GS,
      nPerm = round(1000 / cores) * cores,
      gsSizeLim = c(5, 200),
      ncpus = cores
    )
    cat("Running GSA 3/5")
    gsaRes3 <- runGSA(
      Pval,
      FC,
      geneSetStat = "sum",
      gsc = GS,
      nPerm = round(1000 / cores) * cores,
      gsSizeLim = c(5, 200),
      ncpus = cores
    )
    cat("Running GSA 4/5")
    gsaRes4 <- runGSA(
      Pval,
      FC,
      geneSetStat = "stouffer",
      gsc = GS,
      nPerm = round(1000 / cores) * cores,
      gsSizeLim = c(5, 200),
      ncpus = cores
    )
    cat("Running GSA 5/5")
    gsaRes5 <- runGSA(
      Pval,
      FC,
      geneSetStat = "tailStrength",
      gsc = GS,
      nPerm = round(1000 / cores) * cores,
      gsSizeLim = c(5, 200),
      ncpus = cores
    )
    # No maxmean and fisher: don't support direction. No gsea or
    # wilcoxon, too slow.
    cat("Reorganizing data and prepare for plotting")
    sumTable <-
      GSAsummaryTable(gsaRes2, save = F)  # Needed to extract genesets names and gene numbers
    # Combine results in list
    resList <- list(gsaRes1, gsaRes2, gsaRes3, gsaRes4, gsaRes5)
    resList <-
      setNames(resList,
               c("mean", "median", "sum", "stouffer",
                 "tailStrength"))
    
return(resList)
  }