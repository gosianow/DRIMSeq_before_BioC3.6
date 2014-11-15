##############################################################################
# calculate profile likelihood + adjustements for common dispersion
# dmAdjustedProfileLik, dmSQTLAdjustedProfileLik
##############################################################################
# returns common likelihood = sum of all gene likelihoods



dmSQTLAdjustedProfileLik <- function(gamma0, dgeSQTL, adjust = FALSE, mode = "constrOptim2", epsilon = 1e-05, maxIte = 1000, mcCores=20, verbose = FALSE){
  
  cat("Gamma in optimize:", gamma0, fill = TRUE)
  
  dgeSQTLFit <- dmSQTLFit(dgeSQTL, model = "full", dispersion=gamma0, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose=verbose, mcCores = mcCores)
  
  loglik <- unlist(lapply(dgeSQTLFit$fit, function(g){g$logLik})) 

    loglik <- sum(loglik)
  
  cat("loglik:", loglik, fill = TRUE)
  
  if(!adjust)
    return(loglik)
  
  
  ## Cox-Reid adjustement
  adj <- dmSQTLAdj(gamma0, dgeSQTL = dgeSQTLFit, mcCores=mcCores)

    adj <- sum(adj)
  
  adjloglik <- loglik - adj
  
  cat("adjloglik:", adjloglik, fill = TRUE)
  
  
  return(adjloglik)
  
}


