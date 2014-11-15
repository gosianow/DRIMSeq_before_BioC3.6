
#######################################################
#  group testing
# dmTest, dmSQTLTest
#######################################################




dmSQTLTest <- function(dgeSQTL, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 20){
  
  fit.full <- dgeSQTL$fit

  ## fit null model
  cat("Fitting null model \n")
  print(proc.time())
  fit.null <- dmSQTLFit(dgeSQTL=dgeSQTL, model = "null", dispersion = dgeSQTL$dispersion, mode=mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose, mcCores = mcCores)$fit
  print(proc.time())
  
  cat("Calculating LR \n")
  LRT <- mclapply(dgeSQTL$SNPs$SNP_id, function(snp){
    # snp = "snp_19_502623"
#     if(verbose) 
      cat("testing SNP: ", snp, fill = TRUE)
    
    if(is.null(fit.null[[snp]]) || is.null(fit.full[[snp]])) 
      return(data.frame(LR=NA, df=NA, PValue=NA, LLfull=NA, LLnull=NA))
    
    LLnull <- fit.null[[snp]]$logLik
    
    LLfull <- sum(fit.full[[snp]]$logLik)
    
    LR <-  2*(LLfull - LLnull)
    
    DFnull <- fit.null[[snp]]$df
    DFfull <- sum(fit.full[[snp]]$df)
    
    df <- DFfull - DFnull
    
    pValue <- pchisq(LR, df = df , lower.tail = FALSE)
    
#     gc()
    
    return(data.frame(LR=LR, df=df, PValue=pValue, LLfull=LLfull, LLnull=LLnull))
    
  }, mc.cores=mcCores)
  
print(proc.time())
  save(fit.null, LRT, file="LRT.RData")
  
  LRT <- do.call(rbind, LRT)
  FDR <- p.adjust(LRT$PValue, method="BH")
  o <- order(LRT$PValue)
  
  table <- data.frame(SNP_id = dgeSQTL$SNPs$SNP_id, gene_id = dgeSQTL$SNPs$gene_id, LR=LRT$LR, df=LRT$df, LLfull=LRT$LLfull, LLnull=LRT$LLnull , PValue=LRT$PValue, FDR=FDR, stringsAsFactors = FALSE)[o,]
  
  fit <- dgeSQTL 
  fit$fit.null <- fit.null
  fit$table <- table
  
  return(fit)
  
  
}


