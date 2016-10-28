# Calculate z-scores from two data sets:
# Input:
#           data      data.frame with normalized data. Columns are replicates, rows are genes/fluxes/etc.
#           sampleIdx list of column indices of samples
#           refIdx    list of column indices of reference
# Output:
#           z         z-scores
## 2016-04-19 Eduard Kerkhoven

calcZ<-function(data,sampleIdx,refIdx ){
  samData<-data[,sampleIdx]
  refData<-data[,refIdx]
  samMean<-rowMeans(samData)
  samVar<-apply(samData,1,var)
  refMean<-rowMeans(refData)
  refVar<-apply(refData,1,var)
  z<-(samMean-refMean)/sqrt(samVar+refVar)
}