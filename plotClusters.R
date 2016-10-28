# Plot cluster graphs
# Input:
#           kmeans    Output from kmeans clustering (should contain kmeans$cluster)
#           normData  Matrix containing normalized data (was used as input for kmeans clustering
# Output:
#           plots
## 2016-10-25 Eduard Kerkhoven

plotClusters<-function(k,normdata){
  # Determine number of clusters
  max(k$cluster)
  out <- data.frame(gene=character(),
                   variable=character(), 
                   value=numeric(),
                   cluster=numeric(),
                   stringsAsFactors=FALSE) 
  for (i in 1:max(k$cluster)){ # Loop through all clusters
    dat<-names(k$cluster[k$cluster==i]) # Genes in cluster
    dat<-data.frame(normdata[dat,]) # Extract normalized data for genes in cluster
    dat$gene<-rownames(dat) # Explicitely add gene names
    dat<-melt(dat,id="gene")
    dat$cluster<-i
    out<-rbind(out,dat)
  }
  head(out)
  out<-out
  ggplot(out,aes(x=variable,y=value,group=gene)) + geom_line() + facet_wrap(~cluster,drop=T)
}