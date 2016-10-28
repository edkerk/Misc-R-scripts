filterNAs<-function(x,replicates,samples){
  x[x==0]<-NA
  nas<-is.na(x) # How many NAs are there?
  r<-replicates # Assign number of replicates to r
  naLibs<-matrix(nrow=nrow(x),ncol=samples)
  for (n in 1:samples){
    naLibs[,n]<-rowSums(nas[,(r*n):((r*n)-r+1)])
  }
  out<-as.data.frame(x[!rowSums(naLibs > 1,na.rm=TRUE),])
  return(out)
}
