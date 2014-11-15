##############################################################################
# adjustements to profile likelihood
##############################################################################

dmSQTLAdjCROneGeneManyGroups <- function(y, group, gamma0, piH){  
  # NULL for filtered genes or genes with one exon
  k <- dim(y)[1]
  if(k <= 1) return(NULL)
  
  group <- as.factor(group)
  ngroups <- nlevels(group)
  lgroups <- levels(group)
  
  adjCR = rep(0, ngroups)

  for(gr in 1:ngroups){
    # gr=1
    igroups <- which(group == lgroups[gr])
    a <- dmAdjCROneGeneGroup(y = y[, igroups, drop = FALSE], gamma0, piH = piH[, lgroups[gr]])
    
    if(is.null(a)) return(NULL)
    adjCR[gr] <- a
    
  }
  
   adjCR <- sum(adjCR)
  
  return(adjCR)
  
}

