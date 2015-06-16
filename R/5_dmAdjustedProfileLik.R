##############################################################################
# calculate profile likelihood + adjustements for common dispersion
##############################################################################
# returns common likelihood = sum of all gene likelihoods


dmAdjustedProfileLik <- function(gamma0, dge, adjustDisp = FALSE, modeProp = "constrOptim2", tolProp = 1e-12, verbose = FALSE, BPPARAM = MulticoreParam(workers=1)){
  
  cat("Gamma in optimize:", gamma0, fill = TRUE)
  
  dgeFit <- dmFit(dge, model = "full", dispersion = gamma0, modeProp = modeProp, tolProp = tolProp, verbose = verbose, BPPARAM = BPPARAM)

  logLik <- sum(unlist(lapply(dgeFit$fitFull, function(g){sum(g$logLik)})) )
  
  cat("logLik:", logLik, fill = TRUE)
  
  if(!adjustDisp)
    return(logLik)

  ## Cox-Reid adjustement
  adj <- dmAdj(gamma0, dge = dgeFit, BPPARAM = BPPARAM)

  
  adjLogLik <- logLik - adj
  
  cat("adjLogLik:", adjLogLik, fill = TRUE)
  
  
  return(adjLogLik)
  
}


