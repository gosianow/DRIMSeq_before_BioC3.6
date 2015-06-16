##############################################################################
# calculate common dispersion 
##############################################################################

# adjustDisp = TRUE; intervalDisp = c(0, 1e+3); tolDisp = 1e-01; modeProp = "constrOptim2G"; tolProp = 1e-12; verbose=FALSE; BPPARAM = BPPARAM

dmSQTLEstimateCommonDisp <- function(dgeSQTL, adjustDisp = TRUE, intervalDisp = c(0, 1e+5), tolDisp = 1e-01, modeProp = "constrOptim2G", tolProp = 1e-12, verbose=FALSE, BPPARAM = MulticoreParam(workers=1)){
	
	cat("Estimating common dispersion.. \n")
	### keep only one SNP per gene
	dgeSQTL_subset <- dgeSQTL
	dgeSQTL_subset$genotypes <- lapply(dgeSQTL_subset$genotypes, function(g){ g[1, , drop = FALSE]})
  
  time <- system.time(optimum <- optimize(f = dmSQTLAdjustedProfileLik, interval = intervalDisp,
                  dgeSQTL = dgeSQTL_subset, adjustDisp = adjustDisp, modeProp = modeProp, tolProp = tolProp, verbose = verbose, BPPARAM = BPPARAM,
                  maximum = TRUE, tol = tolDisp) )
  
  cat("Took ", time["elapsed"], " seconds.\n")
  dgeSQTL$commonDispersion <- as.numeric(optimum$maximum)
  cat("** Connom dispersion: ", dgeSQTL$commonDispersion, fill = TRUE)
  return(dgeSQTL)
  
}


