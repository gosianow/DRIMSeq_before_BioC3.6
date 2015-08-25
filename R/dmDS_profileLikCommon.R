##############################################################################
# calculate profile likelihood + adjustements for common dispersion
##############################################################################
# returns common likelihood = sum of likelihoods for all genes

# gamma0 = 38196.6; counts = x@counts; samples = x@samples; disp_adjust = TRUE; prop_mode = "constrOptimG"; prop_tol = 1e-12; verbose = FALSE; BPPARAM = BiocParallel::MulticoreParam(workers = 10)

dmDS_profileLikCommon <- function(gamma0, counts, samples, disp_adjust = TRUE, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  cat("Gamma in optimize:", gamma0, fill = TRUE)
	
  fit_full <- dmDS_fitOneModel(counts = counts, samples = samples, dispersion = gamma0, model = "full", prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
  
  lik <- sum(fit_full@metadata, na.rm = TRUE) ### liks
  
  cat("lik:", lik, fill = TRUE)
  # message("lik:", lik, "\n")
  
  if(!disp_adjust)
    return(lik)

  ## Cox-Reid adjustement for common dispersion
  cat("Calculating adjustement.. \n")
  time <- system.time(adj <- dmDS_adjustmentCommon(gamma0, counts = counts, samples = samples, pi = fit_full, BPPARAM = BPPARAM))
  cat("Took ", time["elapsed"], " seconds.\n")
  
  adjLik <- lik - adj
  
  cat("adjLik:", adjLik, fill = TRUE)
  
  
  return(adjLik)
  
}
