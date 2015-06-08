##############################################################################
# calculate common dispersion 
##############################################################################

dmSQTLEstimateCommonDisp <- function(dgeSQTL, adjustDisp = TRUE, intervalDisp = c(0, 1e+5), tolDisp = 1e-01, modeProp = "constrOptim2G", tolProp = 1e-12, verbose=FALSE, BPPARAM = MulticoreParam(workers=1)){
	
	### keep only one SNP per gene
	dgeSQTL_subset <- dgeSQTL
	dgeSQTL_subset$genotypes <- lapply(dgeSQTL_subset$genotypes, function(g){ g[1, , drop = FALSE]})
  
  optimum <- optimize(f = dmSQTLAdjustedProfileLik, interval = intervalDisp,
                  dgeSQTL = dgeSQTL_subset, adjustDisp = adjustDisp, modeProp = modeProp, tolProp = tolProp, verbose = verbose, BPPARAM = BPPARAM,
                  maximum = TRUE, tol = tolDisp) 
  
  
  dgeSQTL$commonDispersion <- as.numeric(optimum$maximum)
  cat("** Connom dispersion: ", dgeSQTL$commonDispersion, fill = TRUE)
  return(dgeSQTL)
  
}


