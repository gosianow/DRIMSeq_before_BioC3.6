##############################################################################
# adjustements to profile likelihood
##############################################################################


dmSQTLAdj <- function(gamma0, dgeSQTL, mcCores=20){
  
  y <- dgeSQTL$counts
  
  adj <- mclapply(seq(nrow(dgeSQTL$SNPs)), function(snp){  
    # snp = 1
    if(is.null(dgeSQTL$fit[[snp]])) return(NULL)
    
		  NAs <- !is.na(dgeSQTL$genotypes[snp,]) & !is.na(y[[dgeSQTL$SNPs[snp, "gene_id"]]][1, ])
             
		             y.g <- y[[dgeSQTL$SNPs[snp, "gene_id"]]][, NAs]
             
		             group <- dgeSQTL$genotypes[snp, NAs]
								 
	  group <- as.factor(group)
	  ngroups <- nlevels(group)
	  lgroups <- levels(group)
  
	  igroups <- list()
	  for(gr in 1:ngroups){
	    # gr=2
	    igroups[[lgroups[gr]]] <- which(group == lgroups[gr])
    
	  }

		a <- dmAdjCROneGeneManyGroups(y = y.g, ngroups = ngroups, lgroups = lgroups, igroups = igroups, gamma0 = gamma0, piH = dgeSQTL$fit[[snp]]$piH) 
		
    return(a)
    
  }, mc.cores=mcCores)
  
  adj <- unlist(adj)
  adj <- sum(adj[adj != Inf]) ## some genes have adj = Inf so skipp them
	
  return(adj)
  
}


