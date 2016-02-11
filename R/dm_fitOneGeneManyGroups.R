### one gene, many groups
### If something is wrong, return matrix of NA
### returns stats being only likelihoods for different groups (no df)

# y = counts[[g]]; gamma0 = gamma0[g]

dm_fitOneGeneManyGroups <- function(y, ngroups, lgroups, igroups, gamma0, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE){

  k <- nrow(y)
  
  pi <- matrix(NA, nrow = k, ncol = ngroups, dimnames = list(rownames(y), lgroups))
  lik <- rep(NA, ngroups)
  dev <- rep(NA, ngroups)
  names(lik) <- paste0("lik_", lgroups)
  names(dev) <- paste0("dev_", lgroups)
  stats <- c(lik, dev)
  
  if(is.na(gamma0) || k < 2)
  return(list(pi = pi, stats = stats))
  
  
  for(gr in 1:ngroups){
    # gr = 1
    
    fit_gr <- dm_fitOneGeneOneGroup(y = y[, igroups[[gr]], drop = FALSE], gamma0 = gamma0, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)
    
    if(is.na(fit_gr[["stats"]]["lik"])){
      pi <- matrix(NA, nrow = k, ncol = ngroups, dimnames = list(rownames(y), lgroups))
      lik <- rep(NA, ngroups)
      dev <- rep(NA, ngroups)
      names(lik) <- paste0("lik_", lgroups)
      names(dev) <- paste0("dev_", lgroups)
      stats <- c(lik, dev)
      return(list(pi = pi, stats = stats))  
    }
    
    
    pi[,gr] <- fit_gr[["pi"]] ### pi
    lik[gr] <- fit_gr[["stats"]]["lik"] ### likelihood
    dev[gr] <- fit_gr[["stats"]]["dev"] ### deviance
    
  }
  
  return(list(pi = pi, stats = c(lik, dev))) ### pi and stats can have NAs

}


