##############################################################################
# adjustements to profile likelihood
##############################################################################

dm_adjustmentOneGeneManyGroups <- function(y, ngroups, lgroups, igroups, gamma0, pi){  
  # NULL for filtered genes or genes with one exon
  k <- dim(y)[1]
  if(k <= 1) return(NULL)

  adj = rep(0, ngroups)
  
  for(gr in 1:ngroups){
# gr=1
    a <- dm_adjustmentOneGeneOneGroup(y = y[, igroups[[gr]], drop = FALSE], gamma0, pi = pi[, lgroups[gr]])
    
    if(is.null(a)) return(NULL)
			
    adj[gr] <- a
    
  }
  
  adj <- sum(adj)
  
  return(adj)
  
}

