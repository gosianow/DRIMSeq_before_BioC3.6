##############################################################################
# adjustements to profile likelihood
##############################################################################

dmAdjCROneGeneManyGroups <- function(y, ngroups, lgroups, igroups, gamma0, piH){  
  # NULL for filtered genes or genes with one exon
  k <- dim(y)[1]
  if(k <= 1) return(NULL)

  adjCR = rep(0, ngroups)
  
  for(gr in 1:ngroups){
# gr=1
    a <- dmAdjCROneGeneGroup(y = y[, igroups[[gr]], drop = FALSE], gamma0, piH = piH[, lgroups[gr]])
    
    if(is.null(a)) return(NULL)
    adjCR[gr] <- a
    
  }
  
  adjCR <- sum(adjCR)
  
  return(adjCR)
  
}

