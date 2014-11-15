##############################################################################
# adjustements to profile likelihood
##############################################################################

# dge = dgeFit

dmAdj <- function(gamma0, dge, group=NULL, mcCores=20){
  
  y <- dge$counts
  
  if(is.null(group)) group <- dge$samples$group
  group <- as.factor(group)
  ngroups <- nlevels(group)
  lgroups <- levels(group)
  
  igroups <- list()
  for(gr in 1:ngroups){
    # gr=2
    igroups[[lgroups[gr]]] <- which(group == lgroups[gr])
    
  }
  
  adj <- mclapply(seq(length(y)), function(g){  
    # g = 1084

    if(is.null(dge$fit[[g]])) return(NULL)
    
    a <- dmAdjCROneGeneManyGroups(y = y[[g]], ngroups = ngroups, lgroups = lgroups, igroups = igroups, gamma0 = gamma0, piH = dge$fit[[g]]$piH) 
    
    return(a)
    
  }, mc.cores=mcCores)
  
  adj <- unlist(adj)
  adj <- sum(adj[adj != Inf]) ## some genes have adj = Inf so skipp them
  
  return(adj)
  
}

