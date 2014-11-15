##############################################################################
# calculate common dispersion 
# dmEstimateCommonDisp, dmSQTLEstimateCommonDisp
##############################################################################

dmSQTLEstimateCommonDisp <- function(dgeSQTL, adjust = FALSE, subset=Inf, mode = "constrOptim2", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=20, verbose=FALSE){
  
  dgeSQTLorg <- dgeSQTL

  if(subset != Inf){    
    nSNPs <- nrow(dgeSQTL$SNPs)
    
    ## take random SNP
    SNPSubset <- sample.int(nSNPs, min(subset, nSNPs))
  
    dgeSQTL$SNPs <- dgeSQTL$SNPs[SNPSubset, ]
    dgeSQTL$genotypes <- dgeSQTL$genotypes[SNPSubset, ]
    
    
    
  }
  
  out <- optimize(f = dmSQTLAdjustedProfileLik, interval = interval,
                  dgeSQTL = dgeSQTL, adjust = adjust, mode = mode, epsilon = epsilon, maxIte = maxIte, mcCores = mcCores, verbose = verbose,
                  maximum = TRUE, tol = tol) 
  
  
  dgeSQTLorg$commonDispersion <- out$maximum
  cat("** Connom dispersion: ", out$maximum, fill = TRUE)
  return(dgeSQTLorg)
  
}


