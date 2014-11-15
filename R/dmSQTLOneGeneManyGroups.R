

dmSQTLOneGeneManyGroups <- function(y, group, gamma0, mode = "constrOptim2", epsilon = 1e-05, maxIte = 1000, verbose = FALSE){
  # NULL for filtered genes or genes with one exon
  k <- dim(y)[1]
  if(k <= 1) return(NULL)
  
  group <- as.factor(group)
  ngroups <- nlevels(group)
  lgroups <- levels(group)
  
  piH = matrix(0, nrow = k, ncol = ngroups)
  logLik = rep(NA, ngroups)
  df = rep(0, ngroups)

  for(gr in 1:ngroups){
    # gr=2
    igroups <- which(group == lgroups[gr])
    
    fitGr <- dmOneGeneGroup(y = y[, igroups, drop = FALSE], gamma0 = gamma0, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose=verbose)

    if(is.null(fitGr)) return(NULL)
    
    piH[,gr] <- fitGr$piH
    logLik[gr] <- fitGr$logLik
    df[gr] <- fitGr$df
    
  }
    
  colnames(piH) <- names(df) <- lgroups
  logLik <- sum(logLik)

  return(list(piH = piH, gamma0 = gamma0, logLik = logLik, df=df))
  
  
}

