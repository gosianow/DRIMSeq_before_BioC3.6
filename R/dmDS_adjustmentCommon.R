##############################################################################
# adjustements to profile likelihood for common dispersion -> sum
##############################################################################

dmDS_adjustmentCommon <- function(gamma0, counts, samples, pi, BPPARAM = BiocParallel::MulticoreParam(workers=1)){
  
  gene_list <- names(counts)
	
  group <- samples$group
  ngroups <- nlevels(group)
  lgroups <- levels(group)
  
  igroups <- lapply(lgroups, function(gr){which(group == gr)})
  names(igroups) <- lgroups
  
  adj <- BiocParallel::bplapply(gene_list, function(g){  
		# g = gene_list[1]
    
    if(any(is.na(pi[[g]]))) 
			return(NA)

    a <- dm_adjustmentOneGeneManyGroups(y = counts[[g]], ngroups = ngroups, lgroups = lgroups, igroups = igroups, gamma0 = gamma0, pi = pi[[g]]) 
    
    return(a)
    
  }, BPPARAM = BPPARAM)
  
  adj <- unlist(adj)
  adj <- sum(adj[adj != Inf], na.rm = TRUE) ## some genes have adj = Inf so skipp them
  
  return(adj)
  
}

