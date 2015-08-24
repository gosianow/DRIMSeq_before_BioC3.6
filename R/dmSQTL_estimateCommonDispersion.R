##############################################################################
# calculate common dispersion 
##############################################################################
# counts = x@counts; genotypes = x@genotypes; disp_adjust = TRUE; disp_interval = c(0, 1e+5); disp_tol = 1e+01; prop_mode = "constrOptimG"; prop_tol = 1e-12; verbose=FALSE; BPPARAM = BiocParallel::MulticoreParam(workers = 1)

dmSQTL_estimateCommonDispersion <- function(counts, genotypes, disp_adjust = TRUE, disp_interval = c(0, 1e+5), disp_tol = 1e+01, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose=FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
	
	cat("Estimating common dispersion.. \n")
	
	### keep only one SNP per gene
  inds <- 1:length(genotypes)
  names(inds) <- names(genotypes)
	genotypes <- MatrixList(lapply(inds, function(g){ genotypes[[g]][1, , drop = FALSE]}))
	
  time <- system.time(optimum <- optimize(f = dmSQTL_profileLikCommon, interval = disp_interval,
                  counts = counts, genotypes = genotypes, disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM,
                  maximum = TRUE, tol = disp_tol) )
  
	dispersion <- optimum$maximum
		
  cat("Took ", time["elapsed"], " seconds.\n")
  cat("** Connom dispersion: ", dispersion, fill = TRUE)
	
  return(dispersion)
  
}


