
##############################################################################
# calculate profile likelihood + adjustements for common dispersion
# dmAdjustedProfileLik, dmSQTLAdjustedProfileLik
##############################################################################
# returns common likelihood = sum of all gene likelihoods

# gamma0=38196.6; dge <- dge; group=NULL; adjust = TRUE; mode = "constrOptim2G"; epsilon = 1e-05; maxIte = 1000; mcCores=40; verbose = FALSE


dmAdjustedProfileLik <- function(gamma0, dge, group=NULL, adjust = FALSE, mode = "constrOptim2", epsilon = 1e-05, maxIte = 1000, mcCores=20, verbose = FALSE){
  
  cat("Gamma in optimize:", gamma0, fill = TRUE)
  
  dgeFit <- dmFit(dge, group=group, dispersion=gamma0, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose=verbose, mcCores = mcCores)

  loglik <- sum(unlist(lapply(dgeFit$fit, function(g){g$logLik})) )
  
  cat("loglik:", loglik, fill = TRUE)
  
  if(!adjust)
    return(loglik)

  ## Cox-Reid adjustement
  adj <- dmAdj(gamma0, dge = dgeFit, group = group, mcCores=mcCores)

  
  adjloglik <- loglik - adj
  
  cat("adjloglik:", adjloglik, fill = TRUE)
  
  
  return(adjloglik)
  
}


