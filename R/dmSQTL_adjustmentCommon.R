##############################################################################
# adjustements to profile likelihood
##############################################################################


dmSQTL_adjustmentCommon <- function(gamma0, data, fit_full, BPPARAM = MulticoreParam(workers=1)){
  
  gene_list <- names(data@counts)
  
  adj <- bplapply(gene_list, function(g){  
		# g = geneList[1624]; y = data@counts[[g]]; snps = data@genotypes[[g]]
             
		             y = data@counts[[g]]
		             snps = data@genotypes[[g]]
								 
		             adj <- rep(NA, nrow(snps))
             
		             for(i in 1:nrow(snps)){
		               # i = 1
									 
							     if(is.null(fit_full[[g]][[i]]$pi)) 
							 			next
							 
		               NAs <- is.na(snps[i, ]) | is.na(y[1, ])            
		               yg <- y[, !NAs]             
		               group <- snps[i, !NAs]
		               group <- factor(group)
		               ngroups <- nlevels(group)
		               lgroups <- levels(group)
		               nlibs <- length(group)
               
								   igroups <- lapply(lgroups, function(gr){which(group == gr)})
								   names(igroups) <- lgroups

		adj[i] <- dm_adjustmentOneGeneManyGroups(y = yg, ngroups = ngroups, lgroups = lgroups, igroups = igroups, gamma0 = gamma0, pi = fit_full[[g]][[i]]$pi) 
		
	}
    return(adj)
    
  }, BPPARAM = BPPARAM)
  
  adj <- unlist(adj)
  adj <- sum(adj[adj != Inf], na.rm = TRUE) ## some genes have adj = Inf so skipp them
	
  return(adj)
  
}














