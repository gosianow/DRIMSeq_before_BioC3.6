##############################################################################
# calculate profile likelihood + adjustements for common dispersion
# returns common likelihood = sum of likelihoods from all genes
##############################################################################
# gamma0 = 38196.6

dmSQTL_profileLikCommon <- function(gamma0, counts, genotypes, disp_adjust = TRUE, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = MulticoreParam(workers=1)){
  
  cat("Gamma in optimize:", gamma0, fill = TRUE)
  
  fit_full <- dmSQTL_fitOneModel(counts, genotypes, dispersion = gamma0, model = "full", prop_mode = prop_mode, prop_tol = prop_tol, verbose=verbose, BPPARAM = BPPARAM)
  
  lik <- sum(unlist( lapply(fit_full, function(g) g@statistics[, "lik"]) ), na.rm = TRUE) 
  
  cat("lik:", lik, fill = TRUE)
  
  if(!disp_adjust)
    return(lik)
  
  ## Cox-Reid adjustement
  pi <- lapply(fit_full, function(g) g@proportions)
  
  adj <- dmSQTL_adjustmentCommon(gamma0, counts, genotypes, pi = pi , BPPARAM = BPPARAM)
  
  adjLik <- lik - adj
  
  cat("adjLik:", adjLik, fill = TRUE)
  
  
  return(adjLik)
  
}



