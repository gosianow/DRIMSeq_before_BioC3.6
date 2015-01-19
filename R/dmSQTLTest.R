
#######################################################
#  group testing
# dmTest, dmSQTLTest
#######################################################


# dispersion = c("commonDispersion", "tagwiseDispersion")[2]; mode="constrOptim2G"; epsilon = 1e-05; maxIte = 1000; verbose=FALSE; mcCores = 1



dmSQTLTest <- function(dgeSQTL, dispersion = c("commonDispersion", "tagwiseDispersion")[1] , mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 20){
  
  fit.full <- dgeSQTL$fit

  ## fit null model
  cat("Fitting null model \n")

  dgeSQTL.null <- dmSQTLFit(dgeSQTL, model = "null", dispersion = dispersion, mode=mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose, mcCores = mcCores)
  
  fit.null <- dgeSQTL.null$fit
  
  cat("Calculating LR \n")
  LRT <- mclapply(seq(nrow(dgeSQTL$SNPs)), function(snp){
    # snp = 1
			
#     if(verbose) 
      # cat("testing SNP: ", paste0(dgeSQTL$SNPs[snp, c("SNP_id", "gene_id")], collapse = "-"), fill = TRUE)
		
    if(is.null(fit.null[[snp]]) || is.null(fit.full[[snp]])) 
      return(data.frame(LR=NA, df=NA, PValue=NA, LLfull=NA, LLnull=NA))
    
    LLnull <- fit.null[[snp]]$logLik
    
    LLfull <- sum(fit.full[[snp]]$logLik)
    
    LR <-  2*(LLfull - LLnull)
    
    DFnull <- fit.null[[snp]]$df
    DFfull <- sum(fit.full[[snp]]$df)
    
    # df <- DFfull - DFnull
		df <- DFnull * (length(fit.full[[snp]]$df) - 1)
    
    pValue <- pchisq(LR, df = df , lower.tail = FALSE)
    
#     gc()
    
    return(data.frame(LR=LR, df=df, PValue=pValue, LLfull=LLfull, LLnull=LLnull))
    
  }, mc.cores=mcCores)
  
# save(fit.null, LRT, file="LRT.RData")
  
  LRT <- do.call(rbind, LRT)
  FDR <- p.adjust(LRT$PValue, method="BH")
  o <- order(LRT$PValue)
  
  table <- data.frame(SNP_id = dgeSQTL$SNPs$SNP_id, gene_id = dgeSQTL$SNPs$gene_id, LR=LRT$LR, df=LRT$df, LLfull=LRT$LLfull, LLnull=LRT$LLnull , PValue=LRT$PValue, FDR=FDR, stringsAsFactors = FALSE)[o,]
  
  fit <- dgeSQTL 
  fit$fit.null <- fit.null
  fit$table <- table
  
  return(fit)
  
  
}


