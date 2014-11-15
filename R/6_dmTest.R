
#######################################################
#  group testing
# dmTest, dmSQTLTest
#######################################################


# dge = dgeDM; mode="constrOptim2"; epsilon = 1e-05; maxIte = 1000; verbose=FALSE; mcCores=20

dmTest <- function(dge, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 20){
  
  fit.full <- dge$fit
  group <- as.factor(rep(1, length(dge$samples$group)))
  dge$genes$gene_id <- as.factor(dge$genes$gene_id)
  
  ## fit null model
  fit.null <- dmFit(dge=dge, group=group, dispersion=dge$dispersion, mode=mode, epsilon = epsilon, maxIte = maxIte, verbose= verbose, mcCores = mcCores)$fit
  
  LRT <- mclapply(unique(levels(dge$genes$gene_id)), function(g){
    # g = "g1"
    if(verbose) cat("testing gene: ",g, fill = TRUE)
    
    if(is.null(fit.null[[g]]) || is.null(fit.full[[g]])) 
      return(data.frame(LR=NA, df=NA, PValue=NA, LLfull=NA, LLnull=NA))
    
      LLnull <- fit.null[[g]]$logLik

       LLfull <- sum(fit.full[[g]]$logLik)

      LR <-  2*(LLfull - LLnull)
      
      DFnull <- fit.null[[g]]$df
      DFfull <- sum(fit.full[[g]]$df)
      
      df <- DFfull - DFnull
      
      pValue <- pchisq(LR, df = df , lower.tail = FALSE)
      
      return(data.frame(LR=LR, df=df, PValue=pValue, LLfull=LLfull, LLnull=LLnull))
    
  }, mc.cores=mcCores)
  
  
  LRT <- do.call(rbind, LRT)
  FDR <- p.adjust(LRT$PValue, method="BH")
  o <- order(LRT$PValue)
  
  table <- data.frame(GeneID=unique(levels(dge$genes$gene_id)), LR=LRT$LR, df=LRT$df, LLfull=LRT$LLfull, LLnull=LRT$LLnull , PValue=LRT$PValue, FDR=FDR)[o,]
  
  fit <- dge
  fit$fit.null <- fit.null
  fit$table <- table
  
  return(fit)
  
  
}

