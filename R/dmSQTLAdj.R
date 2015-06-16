##############################################################################
# adjustements to profile likelihood
##############################################################################


dmSQTLAdj <- function(gamma0, dgeSQTL, BPPARAM = MulticoreParam(workers=1)){
  
  geneList <- names(dgeSQTL$counts)
  
  adjList <- bplapply(geneList, function(g){  
		# g = geneList[1624]; y = dgeSQTL$counts[[g]]; snps = dgeSQTL$genotypes[[g]]
             
		             y = dgeSQTL$counts[[g]]
		             snps = dgeSQTL$genotypes[[g]]
								 
		             adj <- rep(NA, nrow(snps))
             
		             for(i in 1:nrow(snps)){
		               # i = 1
									 
							     if(is.null(dgeSQTL$fit[[g]][[i]]$piH)) 
							 			next
							 
		               NAs <- is.na(snps[i, ]) | is.na(y[1, ])            
		               yg <- y[, !NAs]             
		               group <- snps[i, !NAs]
		               group <- factor(group)
		               ngroups <- nlevels(group)
		               lgroups <- levels(group)
		               nlibs <- length(group)
               
		               igroups <- list()
		               for(gr in 1:ngroups){
		                 # gr=2
		                 igroups[[lgroups[gr]]] <- which(group == lgroups[gr])
                 
		               }

		adj[i] <- dmAdjCROneGeneManyGroups(y = yg, ngroups = ngroups, lgroups = lgroups, igroups = igroups, gamma0 = gamma0, piH = dgeSQTL$fit[[g]][[i]]$piH) 
		
	}
    return(adj)
    
  }, BPPARAM = BPPARAM)
  
	
		# table(unlist(lapply(1:length(adjList), function(g){
		# class(adjList[[g]][[1]])
		# })))
		#
		# ss <- unlist(lapply(1:length(adjList), function(g){
		# if(class(adjList[[g]][[1]]) == "logical"){
		# 	print(g)
		# print(adjList[[g]][[1]])
		# }
		# }))
	
	
  adj <- unlist(adjList)
  adj <- sum(adj[adj != Inf], na.rm = TRUE) ## some genes have adj = Inf so skipp them
	
  return(adj)
  
}














