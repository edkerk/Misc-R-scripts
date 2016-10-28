# Description:
#     Makes a plot showing the most significant change genesets, as identified
#     from consensus geneset analysis using Piano. For the highest ranking
#     genesets, a barchart is plotted showing the number of up- and down
#     regulated genes. A p-value cutoff is used to signify what genes have
#     significantly changed expression.
#
# Input:      Pval      named vector of gene-level adjusted P-values
#             FC        named vector of gene-level log2 fold-changes
#             GS        gene-set, as loaded usign Piano
#             rankScore cutoff for non-directional consensus score, whether
#                       genesets should be included in barchart.
#             Pcutoff   cutoff for adjusted P-values of differentially expressed
#                       genes. Recommended is 0.05.
#             savePlot  true if plot should be saved as file, false if not
#             title     short title used in plot and filename
# Output:     plot      the plot
#
# 2016-08-22 EJK: Made some adjustments to automatically resize PDF depending
#                 on number of significant GO terms.
# 2016-05-16 Eduard Kerkhoven (eduardk@chalmers.se)



consGSAplot <-
  function(Pval,
           FC,
           GS,
           rankScore,
           Pcutoff,
           savePlot,
           title) {
    require(ggplot2)
    require(piano)
    require(parallel)
    require(snowfall)
    require(tidyr)
    require(scales)
    
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
      gsSizeLim = c(5, 500),
      ncpus = cores
    )
    cat("Running GSA 2/5")
    gsaRes2 <- runGSA(
      Pval,
      FC,
      geneSetStat = "median",
      gsc = GS,
      nPerm = round(1000 / cores) * cores,
      gsSizeLim = c(5, 500),
      ncpus = cores
    )
    cat("Running GSA 3/5")
    gsaRes3 <- runGSA(
      Pval,
      FC,
      geneSetStat = "sum",
      gsc = GS,
      nPerm = round(1000 / cores) * cores,
      gsSizeLim = c(5, 500),
      ncpus = cores
    )
    cat("Running GSA 4/5")
    gsaRes4 <- runGSA(
      Pval,
      FC,
      geneSetStat = "stouffer",
      gsc = GS,
      nPerm = round(1000 / cores) * cores,
      gsSizeLim = c(5, 500),
      ncpus = cores
    )
    cat("Running GSA 5/5")
    gsaRes5 <- runGSA(
      Pval,
      FC,
      geneSetStat = "tailStrength",
      gsc = GS,
      nPerm = round(1000 / cores) * cores,
      gsSizeLim = c(5, 500),
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
    
    # Extract non-directional, distinct up and down genesets.
    non <- consensusScores(resList, class = "non", plot = F)
    up <-
      consensusScores(resList,
                      class = "distinct",
                      direction = "up",
                      plot = F)
    dn <-
      consensusScores(resList,
                      class = "distinct",
                      direction = "down",
                      plot = F)
    non <-
      non$rankMat[non$rankMat[, "ConsScore"] < rankScore + 1, , drop = F]
    dn <-
      dn$rankMat[dn$rankMat[, "ConsScore"] < rankScore + 1, , drop = F]
    up <-
      up$rankMat[up$rankMat[, "ConsScore"] < rankScore + 1, , drop = F]
    dn <- dn[rownames(dn) %in% rownames(non), , drop = F]
    up <- up[rownames(up) %in% rownames(non), , drop = F]
    
    df <- data.frame(Name = c(rownames(up), rownames(dn)))
    df <-
      merge(df, sumTable[, c("Name", "Genes (up)", "Genes (down)",
                             "Genes (tot)")])
    colnames(df) <- c("geneset", "up", "dn", "tot")
    
    df$highup <- 0
    df$highdown <- 0
    df$lowup <- 0
    df$lowdown <- 0
    
    tmp <- gsaRes1$gsc[names(gsaRes1$gsc) %in% df$geneset]
    tmp <- tmp[order(names(tmp))]
    
    # Changed from names(Pval[Pval < Pcutoff]) to rownames(Pval[Pval < Pcutoff,,drop=FALSE]),
    # to allow to work with HUi's proteomics data. Not sure what the difference is.
    for (gset in 1:length(tmp)) {
      df[gset,]$highup <- sum(tmp[[gset]] %in% rownames(Pval[Pval <
                                                            Pcutoff,,drop=FALSE]) &
                                tmp[[gset]] %in% rownames(FC[FC > 0,,drop=FALSE]))  # Highly significant up
      df[gset,]$highdown <- sum(tmp[[gset]] %in% rownames(Pval[Pval <
                                                              Pcutoff,,drop=FALSE]) &
                                  tmp[[gset]] %in% rownames(FC[FC < 0,,drop=FALSE]))  # Highly significant down
      df[gset,]$lowup <- sum(tmp[[gset]] %in% rownames(Pval[Pval >
                                                           Pcutoff,,drop=FALSE]) &
                               tmp[[gset]] %in% rownames(FC[FC > 0,,drop=FALSE]))  # Low significant up
      df[gset,]$lowdown <- sum(tmp[[gset]] %in% rownames(Pval[Pval >
                                                             Pcutoff,,drop=FALSE]) &
                                 tmp[[gset]] %in% rownames(FC[FC < 0,,drop=FALSE]))  # Low significant down
    }
    
    df$geneset <- factor(df$geneset, levels = df[order(df$up / df$dn),
                                                 1])  # Order from mostly up to mostly down.
    
    df2 <- gather(df[, c("geneset", "highup", "highdown", "lowup",
                         "lowdown")], "directAndSignif", "genes", 2:5)
    
    df2$geneset <-
      factor(df2$geneset, levels = df2[order(df$up / df$dn),
                                       1])  # Order from mostly up to mostly down.
    df2$directAndSignif <-
      factor(df2$directAndSignif,
             levels = c("highup",
                        "lowup", "lowdown", "highdown")) # Pointless, order is not maintained by ggplot if stat="identity" is used
    
    # Very inelegant way to order significance and direction of changed genes, instead of order mentioned above
    df2$order <- 0
    df2$order[df2$directAndSignif == "highup"] <- 1
    df2$order[df2$directAndSignif == "lowup"] <- 2
    df2$order[df2$directAndSignif == "lowdown"] <- 3
    df2$order[df2$directAndSignif == "highdown"] <- 4
    df2 <- df2[order(df2$order),]
    
    grht <-
      3 + (dim(df)[1] * 0.5) # Determine height of output graph in cm, 3 cm and additional 0.5 per GO term
    
    ggplot(df2, aes(x = geneset, y = genes, fill = directAndSignif)) +
      geom_bar(position = "fill", stat = "identity") +
      scale_fill_manual(
        values = c("#963836", "#b57372", "#abbdd2", "#4682b4"),
        name = "Direction",
        labels = c(
          paste0("Up (p<", Pcutoff, ")"),
          paste0("Up (p>", Pcutoff, ")"),
          paste0("Down (p>", Pcutoff, ")"),
          paste0("Down (p<", Pcutoff, ")")
        )
      ) 
    
    plot <-
      
      
      ggplot(df2, aes(x = geneset, y = genes, fill = directAndSignif)) +
      geom_bar(position = "fill", stat = "identity") + coord_flip() +
      scale_y_continuous(labels = percent_format(), breaks = c(0,
                                                               0.25, 0.5, 0.75, 1)) +
      scale_fill_manual(
        values = c("#963836", "#b57372", "#abbdd2", "#4682b4"),
        name = "Direction",
        labels = c(
          paste0("Up (p<", Pcutoff, ")"),
          paste0("Up (p>", Pcutoff, ")"),
          paste0("Down (p>", Pcutoff, ")"),
          paste0("Down (p<", Pcutoff, ")")
        )
      ) + labs(x = "GO terms", y = "Genes") +
      theme_set(theme_bw()) + theme(
        text = element_text(size = 7),
        legend.position = "top",
        legend.key.size = unit(7, "points"),
        axis.text = element_text(colour = "black"),
        panel.grid = element_blank(),
        line = element_line(size = 0.25),
        plot.margin = unit(c(1,1, 1, 1), "points"),
        legend.key=element_blank(),
        panel.border=element_blank(),
        axis.ticks.y=element_blank()
      ) + geom_text(
        data = df,
        aes(x = geneset,
            y = 1.04, label = tot),
        inherit.aes = FALSE,
        size = 2
      ) +
      ggtitle(title)
    
    
    if (savePlot == T) {
      plot + ggsave(
        file = paste0("consGSAplot_", title, ".pdf"),
        device = "pdf",
        height = grht,
        units = "cm"
      )
    }
    return(plot)
  }