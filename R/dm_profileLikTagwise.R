##############################################################################
### Tagwise profile likelihood to be optimized
##############################################################################

dm_profileLikTagwise <- function(gamma0, y, ngroups, lgroups, igroups, disp_adjust = TRUE, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE){

  fit <- dm_fitOneGeneManyGroups(y = y, ngroups = ngroups, lgroups = lgroups, igroups = igroups, gamma0 = gamma0, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose) ### return NULL when at least one group has not enough of data 
  
  if(is.null(fit))
    return(NULL)
  
  lik <- sum(fit$lik)
  
  if(!disp_adjust)
    return(lik)
  
  adj <- dm_adjustmentOneGeneManyGroups(y = y, ngroups = ngroups, lgroups = lgroups, igroups = igroups, gamma0 = gamma0, pi = fit$pi) 
  
	if(verbose)
		cat("adj", adj, "\n")
	
  if(adj == Inf)
    return(lik) ### can not return NULL or NA because optim can not handle it
  
  adjLik <- lik - adj
	
#   cat("lik", lik, fill = TRUE)
#   cat("adjLik", adjLik, fill = TRUE)

if(verbose)
	cat("adjLik", adjLik, "\n")
  
  return(adjLik)
  
}

