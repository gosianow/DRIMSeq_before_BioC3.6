#######################################################
#  group testing
#######################################################
# dispersion = c("commonDispersion", "tagwiseDispersion")[2]; modeProp="constrOptim2G"; tolProp = 1e-12; verbose=FALSE

dmSQTLTest <- function(dgeSQTL, dispersion = c("commonDispersion", "tagwiseDispersion")[1] , modeProp="constrOptim2G", tolProp = 1e-12, verbose=FALSE, BPPARAM = MulticoreParam(workers=1)){
  
  # fitFull <- dgeSQTL$fitFull

  ## fit null model
  cat("Fitting null model.. \n")

  time <- system.time(dgeSQTL <- dmSQTLFit(dgeSQTL, model = "null", dispersion = dispersion, modeProp = modeProp, tolProp = tolProp, verbose = verbose, BPPARAM = BPPARAM))
	
  cat("Took ", time["elapsed"], " seconds.\n")
  # fitNull <- dgeSQTL$fitNull
  

  ## calculate LR
  cat("Calculating LR.. \n")
  
  geneList <- names(dgeSQTL$counts)
  
  time <- system.time(LRList <- bplapply(geneList, function(g){
    # g = "ENSG00000188822.6"
    
    fitFull_g <- dgeSQTL$fitFull[[g]]
    fitNull_g <- dgeSQTL$fitNull[[g]]
    
    outTest <- matrix(NA, length(fitFull_g), 6)
    colnames(outTest) <- c("LLfull", "LLnull", "nrGroups", "df", "LR", "pValue")
    
    for(i in 1:length(fitFull_g)){
      # i = 1
      if(is.null(fitFull_g[[i]]) || is.null(fitNull_g[[i]])) 
         next # outTest[i, ] <- NA
      
      outTest[i, "LLfull"] <- sum(fitFull_g[[i]]$logLik)
      outTest[i, "LLnull"] <- fitNull_g[[i]]$logLik
      outTest[i, "nrGroups"] <- length(fitFull_g[[i]]$df)
      outTest[i, "df"] <- fitNull_g[[i]]$df * ( outTest[i, "nrGroups"] - 1 )

    }
    
    NAs <- complete.cases(outTest[, c("LLfull", "LLnull", "nrGroups", "df"), drop = FALSE])
    
    # LR <-  2 * (LLfull - LLnull)
    outTest[NAs, "LR"] <- 2 * (outTest[NAs, "LLfull", drop = FALSE] - outTest[NAs, "LLnull", drop = FALSE])
    outTest[NAs, "pValue"] <- pchisq(outTest[NAs, "LR", drop = FALSE], df = outTest[NAs, "df", drop = FALSE] , lower.tail = FALSE)
    
    return(outTest)
    
  }, BPPARAM = BPPARAM))
	
	cat("Took ", time["elapsed"], " seconds.\n")
	cat("Generating table with results.. \n")
	
	### gene and snp IDs
  IDList <- bplapply(geneList, function(g){
    # g = geneList[1] 
		
		matrix( c(rep(g, length(dgeSQTL$fitFull[[g]])), names(dgeSQTL$fitFull[[g]])), nrow = length(dgeSQTL$fitFull[[g]]), byrow = FALSE, dimnames = list(NULL, c("geneID", "snpID")))
	 
  }, BPPARAM = BPPARAM)
  

	ID <- do.call(rbind, IDList)	
  LR <- data.matrix(do.call(rbind, LRList))
  FDR <- p.adjust(LR[, "pValue"], method="BH")
	
	table <- data.frame(ID, LR, FDR, stringsAsFactors = FALSE)
	
  o <- order(table[, "pValue"])  
  
  dgeSQTL$table <- table[o,]
  
  return(dgeSQTL)
  
  
}


