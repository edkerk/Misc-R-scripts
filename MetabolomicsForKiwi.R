# Generates files for use in Kiwi, based on metabolomics data.
#     sampledata:   Matrix with sample intensities
#     controldata:  Matrix with control intensities
#     prefix: prefix for all files to be generated

MetabolomicsForKiwi<-function(sampledata, controldata, prefix){
  if (!identical(rownames(sampledata),rownames(controldata))){
    stop("Sample data and control data don't have the same rownames, are they derived from the same dataframe?")
    }
  tstats<- apply(cbind(sampledata,controldata), 1, function (x) t.test(x[1:3],x[4:6]))
  tstats<- unlist(lapply(tstats, function(x) x$p.value))
  tstats <-p.adjust(tstats, method = "fdr")
  
  log2fc<- apply(cbind(sampledata,controldata), 1, function (x) log2(mean(x[1:3])/mean(x[4:6])))
  gss<-data.frame(Name=rownames(sampledata))
  gss$`p-value`<-tstats
  rownames(gss[na.omit(gss),])
  
  gls<-gss
  gls$FC<-log2fc
  colnames(gls)<-c("Name","p","FC")
  is.na(gls)
  
  gsc<-cbind(rownames(sampledata),rownames(sampledata))
  
  write.table(gss,file=paste0(prefix,"_GSS.txt"),sep="\t",row.names=F,quote=F)
  write.table(gls,file=paste0(prefix,"_GLS.txt"),sep="\t",row.names=F,quote=F)
  write.table(gsc,file=paste0(prefix,"_GSC.txt"),sep="\t",row.names=F,quote=F)
}

## MANUALLY REMOVE NA VALUES FROM GSS AND GLS FILES (NOT IMPLEMENTED)