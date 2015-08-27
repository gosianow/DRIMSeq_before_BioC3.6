##############################################################################
# calculate profile likelihood + adjustements for common dispersion
# returns common likelihood = sum of likelihoods from all genes
##############################################################################

# gamma0 = 38196.6; counts = x@counts; genotypes = x@genotypes; disp_adjust = TRUE; prop_mode = "constrOptimG"; prop_tol = 1e-12; verbose = FALSE; BPPARAM = BiocParallel::MulticoreParam(workers = 10)


dmSQTL_profileLikCommon <- function(gamma0, counts, genotypes, disp_adjust = TRUE, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  if(verbose) cat("Gamma in optimize:", gamma0, fill = TRUE)
  
  fit_full <- dmSQTL_fitOneModel(counts = counts, genotypes = genotypes, dispersion = gamma0, model = "full", prop_mode = prop_mode, prop_tol = prop_tol, verbose=verbose, BPPARAM = BPPARAM)

  lik <- sum(unlist(lapply(fit_full, function(g) sum(g@metadata, na.rm = TRUE))), na.rm = TRUE)
  
  if(verbose) cat("lik:", lik, fill = TRUE)
  
  if(!disp_adjust)
    return(lik)
  
  ## Cox-Reid adjustement
  if(verbose) cat("* Calculating adjustement.. \n")
  time <- system.time(adj <- dmSQTL_adjustmentCommon(gamma0, counts, genotypes, pi = fit_full, BPPARAM = BPPARAM))
  if(verbose) cat("Took ", time["elapsed"], " seconds.\n")
  
  adjLik <- lik - adj
  
  if(verbose) cat("adjLik:", adjLik, fill = TRUE)
  
  
  return(adjLik)
  
}



