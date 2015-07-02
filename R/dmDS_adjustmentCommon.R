##############################################################################
# adjustements to profile likelihood for common dispersion -> sum
##############################################################################

dmDS_adjustmentCommon <- function(gamma0, data, fit_full, BPPARAM = MulticoreParam(workers=1)){
  
  gene_list <- names(data$counts)
	
  group <- data$samples$group
  ngroups <- nlevels(group)
  lgroups <- levels(group)
  
  igroups <- lapply(lgroups, function(gr){which(group == gr)})
  names(igroups) <- lgroups
  
  adj <- bplapply(gene_list, function(g){  
		
    if(is.null(fit_full[[g]])) 
			return(NULL)

    a <- dm_adjustmentOneGeneManyGroups(y = data$counts[[g]], ngroups = ngroups, lgroups = lgroups, igroups = igroups, gamma0 = gamma0, pi = fit_full[[g]]$pi) 
    
    return(a)
    
  }, BPPARAM = BPPARAM)
  
  adj <- unlist(adj)
  adj <- sum(adj[adj != Inf], na.rm = TRUE) ## some genes have adj = Inf so skipp them
  
  return(adj)
  
}

