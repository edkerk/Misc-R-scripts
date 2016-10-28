# Calculate correlation scores between two conditions. Flux, RNA and protein data are compared.
# Input:
#           geneRxnId   Data.frame with two columns: gene names and rxn IDs
#           flux        Data.frame with four columns of flux data:
#						(1) rxn ID ('rxns')
#						(2) mean reference flux ('refFlux')
#						(3) mean sample flux ('sampleFlux')
#						(4) Z-scores of flux changes. ('ZF')
#           protein     Data.frame with two columns column:
#						(1) genes ('genes')
#						(2) Z-scores of protein changes. ('ZP')
#           RNA         Data.frame with one column:
#						(1) genes ('genes')
#						(2) Z-scores of RNA changes. ('ZR')
# Output:
#           out         data.frame with 6 columns: (1) gene; (2) rxnId; (3) rhoFP-score for
#                       correlation between flux and protein; (4) rhoFR-score for correlation
#                       between flux and RNA; (5) rhoPR-score for correlation between protein
#                       and RNA; and rho-FRP for correlation between flux, protein and RNA.
#                       NaN indicate were values from one of the datasets were missing.
#                       fluxDirectChange indicates that the direction of the flux changed (positive
#                       to negative, or vice versa), which would render the correlation with RNA
#                       and/or protein pointless --> the flux is metabolically regulated.
## 2016-04-27 Eduard Kerkhoven (eduardk@chalmers.se)


# corrFPR<-function(geneRxnId,flux,protein,RNA){
corrFPR<-function(geneRxnId,flux,RNA){
  out<-geneRxnId
  out<-merge(out,flux,by='rxns',sort=F,all.x=T)
  out<-merge(out,RNA,by='genes',sort=F,all.x=T)
#  out<-merge(out,protein,by='genes',sort=F,all.x=T)
  # If reference and sample flux are both negative, then take the inverse Z-score, otherwise keep Z-score
  Idx<-which(out$sampleFlux<0 & out$refFlux<0)
  out$ZF[Idx]<--out$ZF[Idx]
  
  # If the flux changes direction, correlation between flux and protein/RNA doesn't make sense.
#  out$rhoFP[out$sampleFlux*out$refFlux<0]<-'fluxDirectChange' 
  out$rhoFR[out$sampleFlux*out$refFlux<0]<-'fluxDirectChange'
  
  # Calculate rho scores as explained in Bordel et al., 2010
  out$rhoFR<-ifelse(!out$rhoFR%in%'fluxDirectChange',pnorm(out$ZF,0,1)*pnorm(out$ZR,0,1),out$rhoFR)
#  out$rhoFP<-ifelse(!out$rhoFR%in%'fluxDirectChange',pnorm(out$ZF,0,1)*pnorm(out$ZP,0,1),out$rhoFP)
#  out$rhoPR<-pnorm(out$ZR,0,1)*pnorm(out$ZP,0,1)
  
  # If flux and protein/RNA are both down, then use absolute Z-scores to calculate rho-scores
  out$rhoFR<-ifelse(out$ZF<0 & out$ZR < 0 & !out$rhoFR%in%'fluxDirectChange', pnorm(abs(out$ZF),0,1)*pnorm(abs(out$ZR),0,1),out$rhoFR)
#  out$rhoFP<-ifelse(out$ZF<0 & out$ZP < 0 & !out$rhoFP%in%'fluxDirectChange', pnorm(abs(out$ZF),0,1)*pnorm(abs(out$ZP),0,1),out$rhoFP)
#  out$rhoPR<-ifelse(out$ZR < 0 & out$ZP < 0,pnorm(abs(out$ZR),0,1)*pnorm(abs(out$ZP),0,1),out$rhoPR)
  
  # Calculate final rhoFPR score by multiplying rhoFP and rhoFR scores.
 # out$rhoFPR<-suppressWarnings(t(as.numeric(out$rhoFP)*t(as.numeric(out$rhoFR)))) # Suppress warning that NAs are introduced due to multiplication of NaN values
  return(out)
}