##############################################################################
# calculate common dispersion 
##############################################################################

dmEstimateCommonDisp <- function(dge, adjustDisp = TRUE, intervalDisp = c(0, 1e+5), tolDisp = 1e-00, modeProp = "constrOptim2G", tolProp = 1e-12, verbose=FALSE, BPPARAM = MulticoreParam(workers=1)){
	
  cat("Estimating common dispersion.. \n")

  time <- system.time(optimum <- optimize(f = dmAdjustedProfileLik, interval = intervalDisp,
                  dge = dge, adjustDisp = adjustDisp, modeProp = modeProp, tolProp = tolProp, verbose = verbose, BPPARAM = BPPARAM,
                  maximum = TRUE, tol = tolDisp) )
  
  cat("Took ", time["elapsed"], " seconds.\n")
  dge$commonDispersion <- optimum$maximum
  cat("** Connom dispersion: ", optimum$maximum, fill = TRUE)
  return(dge)
  
}

