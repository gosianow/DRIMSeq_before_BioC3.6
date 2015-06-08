#######################################################
#  group testing
#######################################################


dmSQTLTest <- function(dgeSQTL, dispersion = c("commonDispersion", "tagwiseDispersion")[1] , modeProp="constrOptim2G", tolProp = 1e-12, verbose=FALSE, BPPARAM = MulticoreParam(workers=1)){
  
  fitFull <- dgeSQTL$fitFull

  ## fit null model
  cat("Fitting null model.. \n")

  dgeSQTL <- dmSQTLFit(dgeSQTL, model = "null", dispersion = dispersion, modeProp = modeProp, tolProp = tolProp, verbose = verbose, BPPARAM = BPPARAM)
  
  fitNull <- dgeSQTL$fitNull
  
  
  ## calculate LR
  cat("Calculating LR.. \n")
  
  geneList <- names(dgeSQTL$counts)
  
  LRList <- bplapply(geneList, function(g){
    # g = geneList[1]
    
    fitFull_g <- fitFull[[g]]
    fitNull_g <- fitNull[[g]]
    
    outTest <- matrix(0, length(fitFull_g), 6)
    colnames(outTest) <- c("LLfull", "LLnull", "nrGroups", "df", "LR", "pValue")
    
    for(i in 1:length(fitFull_g)){
      # i = 1
      if(is.null(fitFull_g[[i]]) || is.null(fitNull_g[[i]])) 
        outTest[i, ] <- NA
      
      outTest[i, "LLfull"] <- sum(fitFull_g[[i]]$logLik)
      outTest[i, "LLnull"] <- fitNull_g[[i]]$logLik
      outTest[i, "nrGroups"] <- length(fitFull_g[[i]]$df)
      outTest[i, "df"] <- fitNull_g[[i]]$df * ( outTest[i, "nrGroups"] - 1 )

    }
    
    NAs <- complete.cases(outTest)
    
    # LR <-  2*(LLfull - LLnull)
    outTest[NAs, "LR"] <- 2 * (outTest[NAs, "LLfull"] - outTest[NAs, "LLnull"])
    outTest[NAs, "pValue"] <- pchisq(outTest[NAs, "LR"], df = outTest[NAs, "df"] , lower.tail = FALSE)
    
    LR <- data.frame(geneID = g, snpID = names(fitFull_g), outTest, stringsAsFactors = FALSE)
    
    return(LR)
    
  }, BPPARAM = BPPARAM)
  
  table <- do.call(rbind, LRList)
  table$FDR <- p.adjust(table$pValue, method="BH")
  o <- order(table$pValue)
  
  table <- table[o,]
  
  dgeSQTL$table <- table
  
  return(dgeSQTL)
  
  
}


