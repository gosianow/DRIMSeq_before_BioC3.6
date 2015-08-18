##############################################################################
# calculate common dispersion 
##############################################################################

dmDS_estimateCommonDispersion <- function(counts, samples, disp_adjust = TRUE, disp_interval = c(0, 1e+5), disp_tol = 1e-01, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose=FALSE, BPPARAM = BiocParallel::MulticoreParam(workers=1)){
	
  cat("Estimating common dispersion.. \n")

  time <- system.time(optimum <- optimize(f = dmDS_profileLikCommon, interval = disp_interval,
                  counts = counts, samples = samples, disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM,
                  maximum = TRUE, tol = disp_tol) )
  
									dispersion <- optimum$maximum
									
  cat("Took ", time["elapsed"], " seconds.\n")
  cat("** Connom dispersion: ", dispersion, fill = TRUE)
	
  return(dispersion)
  
}

