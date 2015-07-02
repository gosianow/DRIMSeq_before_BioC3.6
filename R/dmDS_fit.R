##############################################################################
# multiple group fitting 
##############################################################################


dmDS_fit <- function(data, dispersion, prop_mode = c("constrOptim2", "constrOptim2G", "FisherScoring")[2], prop_tol = 1e-12, verbose = FALSE, BPPARAM = MulticoreParam(workers=1)){
  
		fit <- list()
	
	fit$fit_full <- dmDS_fitOneModel(data = data, dispersion = dispersion, model = "full", prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)  
  
	fit$fit_null <- dmDS_fitOneModel(data = data, dispersion = dispersion, model = "null", prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
  
	return(fit)
	
}


