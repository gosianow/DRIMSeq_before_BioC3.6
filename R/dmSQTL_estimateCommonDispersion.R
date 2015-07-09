##############################################################################
# calculate common dispersion 
##############################################################################


dmSQTL_estimateCommonDispersion <- function(data, disp_adjust = TRUE, disp_interval = c(0, 1e+5), disp_tol = 1e-01, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose=FALSE, BPPARAM = MulticoreParam(workers=1)){
	
	cat("Estimating common dispersion.. \n")
	
	### keep only one SNP per gene
	genotypes <- lapply(data@genotypes, function(g){ g[1, , drop = FALSE]})  
	data_subset <- new("dmSQTLdata", counts = data@counts, genotypes = genotypes, samples = data@samples)
	
  time <- system.time(optimum <- optimize(f = dmSQTL_profileLikCommon, interval = disp_interval,
                  data = data_subset, disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM,
                  maximum = TRUE, tol = disp_tol) )
  
	  dispersion <- optimum$maximum
		
  cat("Took ", time["elapsed"], " seconds.\n")
  cat("** Connom dispersion: ", dispersion, fill = TRUE)
	
  return(dispersion)
  
}


