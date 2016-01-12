##############################################################################
# calculate common dispersion 
##############################################################################

# counts = x@counts; genotypes = x@genotypes; disp_adjust = TRUE; disp_interval = c(0, 1e+5); disp_tol = 1e+01; prop_mode = "constrOptimG"; prop_tol = 1e-12; verbose = TRUE; BPPARAM = BiocParallel::MulticoreParam(workers = 10)

dmSQTL_estimateCommonDispersion <- function(counts, genotypes, disp_adjust = TRUE, disp_interval = c(0, 1e+5), disp_tol = 1e+01, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  if(verbose) cat("* Estimating common dispersion.. \n")
  
  time <- system.time(optimum <- optimize(f = dmSQTL_profileLikCommon, interval = disp_interval,
                  counts = counts, genotypes = genotypes, disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, verbose = max(0, verbose-1), BPPARAM = BPPARAM,
                  maximum = TRUE, tol = disp_tol) )
  
  dispersion <- optimum$maximum
    
  if(verbose) cat("Took ", time["elapsed"], " seconds.\n")
  if(verbose) cat("*** Connom dispersion: ", dispersion, fill = TRUE)
  
  return(dispersion)
  
}


