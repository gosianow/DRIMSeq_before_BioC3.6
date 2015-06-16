#######################################################
#  group testing
#######################################################
# dge = dgeDM; dispersion = c("commonDispersion", "tagwiseDispersion")[2]; modeProp="constrOptim2G"; tolProp = 1e-12; verbose=FALSE

dmTest <- function(dge, dispersion = c("commonDispersion", "tagwiseDispersion")[1] , modeProp="constrOptim2G", tolProp = 1e-12, verbose=FALSE, BPPARAM = MulticoreParam(workers=1)){
  
  # fitFull <- dge$fitFull
	
  ## fit null model
	cat("Fitting null model.. \n")
	
   time <- system.time(dge <- dmFit(dge = dge, model = "null", dispersion = dispersion, modeProp = modeProp, tolProp = tolProp, verbose = verbose, BPPARAM = BPPARAM))
	 
	 cat("Took ", time["elapsed"], " seconds.\n")
	# fitNull <- dge$fitNull
	
  ## calculate LR
  cat("Calculating LR.. \n")
  
  geneList <- names(dge$counts)
	
  time <- system.time(LRList <- bplapply(geneList, function(g){
    # g = geneList[1]
    if(verbose) cat("testing gene: ", g, fill = TRUE)
    
    if(is.null(dge$fitNull[[g]]) || is.null(dge$fitFull[[g]])) 
    return(rep(NA, 6))
		
      LLnull <- dge$fitNull[[g]]$logLik

      LLfull <- sum(dge$fitFull[[g]]$logLik)

      LR <-  2*(LLfull - LLnull)
			
      nrGroups <- length(dge$fitFull[[g]]$df)
			
      # df <- DFfull - DFnull
      df <- dge$fitNull[[g]]$df * (nrGroups - 1) # (k-1) * nr of groups
			
      pValue <- pchisq(LR, df = df , lower.tail = FALSE)
      
    return(c(LLfull, LLnull, nrGroups, df, LR, pValue))
  }, BPPARAM = BPPARAM))
  
	cat("Took ", time["elapsed"], " seconds.\n")
	cat("Generating table with results.. \n")
	
  LR <- do.call(rbind, LRList)
	colnames(LR) <- c("LLfull", "LLnull", "nrGroups", "df", "LR", "pValue")
  FDR <-  p.adjust(LR[, "pValue"], method="BH")
  
	table <- data.frame(geneID = geneList, LR, FDR, stringsAsFactors = FALSE)
	
	o <- order(table[, "pValue"])
	
  dge$table <- table[o,]
  
  return(dge)
  
  
}

