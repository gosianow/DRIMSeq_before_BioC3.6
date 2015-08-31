##############################################################################
### Tagwise profile likelihood to be optimized
##############################################################################

# gamma0 = splineDisp[i]; y = counts[[g]]

dm_profileLikTagwise <- function(gamma0, y, ngroups, lgroups, igroups, disp_adjust = TRUE, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE){

  fit <- dm_fitOneGeneManyGroups(y = y, ngroups = ngroups, lgroups = lgroups, igroups = igroups, gamma0 = gamma0, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose) ### return NA when at least one group has not enough of data 
  
  if(sum(!is.na(fit$stats)) < 2)
  return(NA)
  
  lik <- sum(fit$stats, na.rm = TRUE)
  
  if(!disp_adjust)
    return(lik)
  
  adj <- dm_adjustmentOneGeneManyGroups(y = y, ngroups = ngroups, lgroups = lgroups, igroups = igroups, gamma0 = gamma0, pi = fit$pi) 
  
# cat("adj", adj, "\n")
	
  adjLik <- lik - adj
	
#   cat("lik", lik, fill = TRUE)
#   cat("adjLik", adjLik, fill = TRUE)

# if(verbose) cat("adjLik", adjLik, "\n")
  
  return(adjLik)
  
}

