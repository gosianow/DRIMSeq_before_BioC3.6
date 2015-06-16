##############################################################################
# calculate profile likelihood + adjustements for common dispersion
# returns common likelihood = sum of likelihoods from all genes
##############################################################################

# adjustDisp = TRUE; modeProp = "constrOptim2G"; tolProp = 1e-12; verbose = FALSE; 

# model = "full"; dispersion = gamma0

dmSQTLAdjustedProfileLik <- function(gamma0, dgeSQTL, adjustDisp = TRUE, modeProp = "constrOptim2G", tolProp = 1e-12, verbose = FALSE, BPPARAM = MulticoreParam(workers=1)){
  
  cat("Gamma in optimize:", gamma0, fill = TRUE)
  
  dgeSQTLFit <- dmSQTLFit(dgeSQTL, model = "full", dispersion = gamma0, modeProp = modeProp, tolProp = tolProp, verbose=verbose, BPPARAM = BPPARAM)
  
  logLik <- sum(unlist( lapply(dgeSQTLFit$fitFull, function(g){lapply(g, function(s){sum(s$logLik)})}) )) 
  
  cat("logLik:", logLik, fill = TRUE)
  
  if(!adjustDisp)
    return(logLik)
  
  ## Cox-Reid adjustement
  adj <- dmSQTLAdj(gamma0, dgeSQTL = dgeSQTLFit, BPPARAM = BPPARAM)
  
  adjLogLik <- logLik - adj
  
  cat("adjLogLik:", adjLogLik, fill = TRUE)
  
  
  return(adjLogLik)
  
}



