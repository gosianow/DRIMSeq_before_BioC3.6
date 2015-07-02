##############################################################################
# calculate profile likelihood + adjustements for common dispersion
##############################################################################
# returns common likelihood = sum of likelihoods for all genes


dmDS_profileLikCommon <- function(gamma0, data, disp_adjust = FALSE, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = MulticoreParam(workers=1)){
  
  cat("Gamma in optimize:", gamma0, fill = TRUE)
	
  fit_full <- dmDS_fitOneModel(data, dispersion = gamma0, model = "full", prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)

  lik <- sum(unlist(lapply(fit_full, function(g){sum(g$lik)})) )
  
  cat("lik:", lik, fill = TRUE)
  
  if(!disp_adjust)
    return(lik)

  ## Cox-Reid adjustement for common dispersion
  adj <- dmDS_adjustmentCommon(gamma0, data = data, fit_full = fit_full, BPPARAM = BPPARAM)

  
  adjLik <- lik - adj
  
  cat("adjLik:", adjLik, fill = TRUE)
  
  
  return(adjLik)
  
}


