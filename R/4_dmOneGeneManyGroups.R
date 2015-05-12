

### -- one gene, many groups

# y <- y[["ENSG00000135778"]]; mode = "constrOptim2"; epsilon = 1e-05; maxIte = 1000; verbose=FALSE

dmOneGeneManyGroups <- function(y, ngroups, lgroups, igroups, gamma0, mode = "constrOptim2", epsilon = 1e-05, maxIte = 1000, verbose = FALSE){
  # NULL for filtered genes or genes with one exon
  k <- dim(y)[1]
  if(k <= 1) return(NULL)
  
  piH = matrix(0, nrow = k, ncol = ngroups)
  logLik = rep(NA, ngroups)
  df = rep(0, ngroups)
  
  for(gr in 1:ngroups){
    # gr=2
    # cat(gr, fill = TRUE)
    
    fitGr <- dmOneGeneGroup(y = y[, igroups[[gr]], drop = FALSE], gamma0 = gamma0, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose=FALSE)
    
    if(is.null(fitGr)) return(NULL)
    
    piH[,gr] <- fitGr$piH
    logLik[gr] <- fitGr$logLik
    df[gr] <- fitGr$df
    
  }
  
  colnames(piH) <- names(df) <- lgroups
  logLik <- sum(logLik)
  
	if(verbose){
		cat("gamma0", gamma0, "\n")
		cat("logLik", logLik, "\n")
	}
	
  return(list(piH = piH, gamma0 = gamma0, logLik = logLik, df=df))
  
}


