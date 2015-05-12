##############################################################################
# calculate tagwise dispersions 
# dmEstimateTagwiseDisp, dmSQTLEstimateTagwiseDisp
##############################################################################

### Tagwise likelihood to be optimized
dmAdjustedProfileLikTG <- function(gamma0, y, ngroups, lgroups, igroups, adjust = FALSE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose = FALSE){
  
#   cat("gamma0",gamma0, fill = TRUE)
  
  f <- dmOneGeneManyGroups(y=y, ngroups=ngroups, lgroups=lgroups, igroups=igroups, gamma0=gamma0, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose) ### return NULL when at least one group has not enought of data 
  
  if(is.null(f))
    return(NULL)
  
  loglik <- f$logLik
  
  if(!adjust)
    return(loglik)
  
  adj <- dmAdjCROneGeneManyGroups(y = y , ngroups = ngroups, lgroups = lgroups, igroups = igroups, gamma0 = gamma0, piH = f$piH) 
  
	if(verbose)
		cat("adj", adj, "\n")
	
  if(adj == Inf)
    return(loglik) ### can not return NULL or NA because optim can not handle it
  
  adjloglik <- loglik - adj
	
#   cat("loglik",loglik, fill = TRUE)
#   cat("adjloglik",adjloglik, fill = TRUE)

if(verbose)
	cat("adjloglik", adjloglik, "\n")
  
  return(adjloglik)
  
}

