### filterLowReads
#
#   This function can filter out low RNAseq reads. For each set of replicates, it check whether each read is over the
#   cpmCutoff threshold (recommended value: 1). It then requires that a minimum of replCutoff (recommended value: 2)
#   replicates have NOT reached this cutoff. For a gene to be filtered out, more than condCutoff (recommended value: 1)
#   conditions should not reached the replCutoff.
#
#   Example:
#   Condition 1       Condition 2       Condition 3
#   #1    #2    #3    #1    #2    #3    #1    #2    #3    #4
#   0     1     12    621   578   547   123   3     0     0
#   12    23    12    0     0     56    145   576   0     0
#
#   First gene:
#   If these are CPM values, 1#1, 3#3 and 3#4 don't reach the cpmCutoff of 1.
#   All conditions reach replCutoff of 2, because for each condition there are at least 2 replicates with enough counts.
#   The first gene is not filtered out, because more than 1 condition has enough replicates with enough counts.
#
#   Second gene:
#   Now 2#1 and 2#2 don't reach the cpmCutoff, and also condition 2 completely fails, because there are less than 2
#   replicates with enough counts. However, the gene is still not filtered out, because 1 condition is allowed to not reach
#   the cutoff.
#
#   2016-12-08    Eduard Kerkhoven    eduardk@chalmers.se

filterLowReads<-function(x,repl,cpmCutoff,replCutoff,condCutoff){
  # cpmCutoff is the minimum CPM (count per million) (recommended: 1)
  # replCutoff is the minimum number of replicates that reach cpmCutoff (recommended: 2)
  # condCutoff is the maximum number of conditions that don't reach replCutoff (recommended: 1)
  require(edgeR)
  cpms<-cpm(x) # Convert to cpm. Requires edgeR
  cols<-matrix(nrow=length(repl),ncol=2) # Make a structure with start and end column for each set of replicates
  cols[1,]<-c(1,repl[1])
  for (i in 2:length(repl)){ # Populate column ID structure. Makes it easier to apply cpm cutoff per set of replicates
    cols[i,1]<-cols[i-1,2]+1
    cols[i,2]<-cols[i,1]+repl[i]-1
  }
  lowLibs<-matrix(nrow=dim(x)[1],ncol=length(repl)) # Make a structure were low reads are indicated
  for (i in 1:length(repl)){
    lowLibs[,i]<-rowSums(cpms[,c(cols[i,1]:cols[i,2])]<cpmCutoff)>replCutoff
  }
  out<-x[!rowSums(lowLibs)>replCutoff,]
}