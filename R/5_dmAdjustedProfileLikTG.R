##############################################################################
### Tagwise likelihood to be optimized
##############################################################################

dmAdjustedProfileLikTG <- function(gamma0, y, ngroups, lgroups, igroups, adjustDisp = TRUE, modeProp = "constrOptim2G", tolProp = 1e-12, verbose = FALSE){
  
#   cat("gamma0", gamma0, fill = TRUE)
  
  f <- dmOneGeneManyGroups(y = y, ngroups = ngroups, lgroups = lgroups, igroups = igroups, gamma0 = gamma0, modeProp = modeProp, tolProp = tolProp, verbose = verbose) ### return NULL when at least one group has not enough of data 
  
  if(is.null(f))
    return(NULL)
  
  logLik <- sum(f$logLik)
  
  if(!adjustDisp)
    return(logLik)
  
  adj <- dmAdjCROneGeneManyGroups(y = y, ngroups = ngroups, lgroups = lgroups, igroups = igroups, gamma0 = gamma0, piH = f$piH) 
  
	if(verbose)
		cat("adj", adj, "\n")
	
  if(adj == Inf)
    return(logLik) ### can not return NULL or NA because optim can not handle it
  
  adjLogLik <- logLik - adj
	
#   cat("loglik",loglik, fill = TRUE)
#   cat("adjloglik",adjloglik, fill = TRUE)

if(verbose)
	cat("adjLogLik", adjLogLik, "\n")
  
  return(adjLogLik)
  
}

