##############################################################################
# calculate common dispersion 
# dmEstimateCommonDisp, dmSQTLEstimateCommonDisp
##############################################################################

# group=NULL; adjust = FALSE; mode = "constrOptim2G"; epsilon = 1e-05; maxIte = 1000; interval = c(0, 1e+5); tol = 1e-00; mcCores=20; verbose=FALSE


dmEstimateCommonDisp <- function(dge, group=NULL, adjust = FALSE, subset=Inf, mode = "constrOptim2", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=20, verbose=FALSE){
  
  if(is.null(group)) group <- dge$samples$group

  dgeOrg <- dge
  
  if(subset != Inf){
    
    allGenes <- unique(as.character(dge$genes$gene_id))
    nGenes <- length(allGenes)
    
    ## take random genes
    genesSubset <- sample.int(nGenes, min(subset, nGenes))
    
    dge <- dge[dge$genes$gene_id %in% allGenes[genesSubset], ]
    dge$genes$gene_id <- as.factor(as.character(dge$genes$gene_id)) # otherwise number of levels stays as before subsetting      
    
  }

  out <- optimize(f = dmAdjustedProfileLik, interval = interval,
                  dge = dge, group = group, adjust = adjust, mode = mode, epsilon = epsilon, maxIte = maxIte, mcCores = mcCores, verbose = verbose,
                  maximum = TRUE, tol = tol) 
  
  
  dgeOrg$commonDispersion <- out$maximum
  cat("** Connom dispersion: ", out$maximum, fill = TRUE)
  return(dgeOrg)
  
}

