##############################################################################
# adjustements to profile likelihood
##############################################################################

# dge = dgeFit

dmAdj <- function(gamma0, dge, BPPARAM = MulticoreParam(workers=1)){
  
  geneList <- names(dge$counts)
	
  group <- factor(dge$samples$group)
  ngroups <- nlevels(group)
  lgroups <- levels(group)
  
  igroups <- list()
  for(gr in 1:ngroups){
    # gr=2
    igroups[[lgroups[gr]]] <- which(group == lgroups[gr])
    
  }
  
  adj <- bplapply(geneList, function(g){  
    # g = 

    if(is.null(dge$fitFull[[g]])) 
			return(NULL)

    a <- dmAdjCROneGeneManyGroups(y = dge$counts[[g]], ngroups = ngroups, lgroups = lgroups, igroups = igroups, gamma0 = gamma0, piH = dge$fitFull[[g]]$piH) 
    
    return(a)
    
  }, BPPARAM = BPPARAM)
  
  adj <- unlist(adj)
  adj <- sum(adj[adj != Inf]) ## some genes have adj = Inf so skipp them
  
  return(adj)
  
}

