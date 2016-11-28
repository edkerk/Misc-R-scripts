# Plot cluster graphs
# Input:
#           clusters  Cluster information per gene (from kmeans: fit$cluster; from hclust: after cutree)
#           normData  Matrix containing normalized data (was used as input for kmeans clustering
#           centroid  Either 'centroid', 'genes' or 'both'
# Output:
#           plots
## 2016-10-25 Eduard Kerkhoven

plotClusters <- function(clusters, normdata, centroid) {
  # Determine number of clusters
    max(clusters)
    out <- data.frame(
      gene = character(),
      variable = character(),
      value = numeric(),
      cluster = numeric(),
      stringsAsFactors = FALSE
    )
    centrOut <- data.frame(
      value=numeric(),
      cluster=numeric(),
      variable=character(),
      stringsAsFactors=FALSE
    )
    for (i in 1:max(clusters)) {
      # Loop through all clusters
      dat <- names(clusters[clusters == i]) # Genes in cluster
      dat <-
        data.frame(normdata[dat, ]) # Extract normalized data for genes in cluster
      cent<-data.frame(value=colMeans(dat,na.rm=T))
      cent$cluster<-i
      cent$variable<-rownames(cent)
      centrOut<-rbind(centrOut,cent)
      dat$gene <- rownames(dat) # Explicitely add gene names
      dat <- melt(dat, id = "gene")
      dat$cluster <- i
      out <- rbind(out, dat)
    }
    head(out)
    out <- out
    switch(centroid,
           both=ggplot(out, aes(x = variable, y = value, group = gene)) + geom_line() + geom_line(data=centrOut, aes(x=variable,y=value,group=1), colour='red',size=1)  + facet_wrap( ~cluster, drop = T),
           centroid=ggplot() + geom_line(data=centrOut, aes(x=variable,y=value,group=1), colour='red',size=1)  + facet_wrap( ~cluster, drop = T),
           genes=ggplot(out, aes(x = variable, y = value, group = gene)) + geom_line() + facet_wrap( ~cluster, drop = T)
    )
}