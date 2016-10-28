# Generate a data.frame comparing two data sets:
# Input:
#           x         Matrix with all the data, first columns are control,
#                     remaining columns are samples
#           repl      Number of replicates, should be equal for both data sets
# Output:
#           out$FC    Log 2 fold-change
#           out$adjP  Adjusted p-value (FDR corrected)
#           out$z     Z-score
## 2015-12-22 Eduard Kerkhoven

tTestTable<-function(x,repl){ # is matrix with all data, repl is number of replicates
  fc<-log2(rowMeans(x[,(repl+1):(2*repl)])/rowMeans(x[,1:repl]))
  pval <- matrix(apply(x,1,function(x){t.test(x[1:repl],x[(repl+1):(repl*2)],paired=F)}$p.value))
  adjpval <- p.adjust(pval,"BY")
  vars1<-apply(x[,1:repl],1,var) # Control
  vars2<-apply(x[,(repl+1):(2*repl)],1,var) # Sample
  z<-(rowMeans(x[,(repl+1):(2*repl)])-rowMeans(x[,1:repl]))/sqrt(vars1+vars2)
  df<-data.frame(log2FC=fc)
  df$p<-pval
  df$adjP<-adjpval
  df$z<-z
  row.names(df)<-row.names(x)
  return(df)
}