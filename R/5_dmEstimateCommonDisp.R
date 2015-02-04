##############################################################################
# calculate common dispersion 
# dmEstimateCommonDisp, dmSQTLEstimateCommonDisp
##############################################################################

# group=NULL; adjust = FALSE; mode = "constrOptim2G"; epsilon = 1e-05; maxIte = 1000; interval = c(0, 1e+5); tol = 1e-00; mcCores=20; verbose=FALSE


dmEstimateCommonDisp <- function(dge, group=NULL, adjust = TRUE, subset=Inf, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=20, verbose=FALSE){
  
  if(is.null(group)) group <- dge$samples$group

	  y <- dge$counts
	  genes <- names(y)
	  ngenes <- length(y)

	  ### calculate mean expression of genes 
	  meanExpr <- unlist(mclapply(seq(ngenes), function(g){ sum(y[[g]]) / ncol(y[[g]]) },  mc.cores=mcCores))  
	  names(meanExpr) <- genes
	  dge$meanExpr <- meanExpr

  dgeOrg <- dge
  
  ### use genes with high gene expression 
  if(subset != Inf){
    
		if(subset >= 1){		    
# 		  ## take random genes
# 		  allGenes <- unique(as.character(dge$genes$gene_id))
# 		  nGenes <- length(allGenes)
# 		  genesSubset <- sample.int(nGenes, min(subset, nGenes))
# 		  dge <- dge[dge$genes$gene_id %in% allGenes[genesSubset], ]     
      
		  minExpr <- subset
		  cat("* Min expression: ", minExpr, "\n")
		  dge <- dge[dge$meanExpr > minExpr, ]
      
	}else{
		minExpr <- quantile(dge$meanExpr, probs = 1 - subset)
		cat("* Min expression: ", minExpr, "\n")
		dge <- dge[dge$meanExpr > minExpr, ]
	}

    dge$genes$gene_id <- factor(as.character(dge$genes$gene_id)) # otherwise number of levels stays as before subsetting      
    
  }

  out <- optimize(f = dmAdjustedProfileLik, interval = interval,
                  dge = dge, group = group, adjust = adjust, mode = mode, epsilon = epsilon, maxIte = maxIte, mcCores = mcCores, verbose = verbose,
                  maximum = TRUE, tol = tol) 
  
  
  dgeOrg$commonDispersion <- out$maximum
  cat("** Connom dispersion: ", out$maximum, fill = TRUE)
  return(dgeOrg)
  
}

