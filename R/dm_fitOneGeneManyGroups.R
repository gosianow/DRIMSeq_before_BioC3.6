### one gene, many groups


dm_fitOneGeneManyGroups <- function(y, ngroups, lgroups, igroups, gamma0, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE){
  
  # If something is wrong, return matrix of NA
  k <- nrow(y)
  
  pi <- matrix(NA, nrow = k, ncol = ngroups, dimnames = list(rownames(y), lgroups))
  stats <- rep(NA, 2)
  names(stats) <- c("lik", "df")
 
  if(is.na(gamma0))
  return(list(pi = pi, stats = stats))
  
  if(k <= 1) 
  return(list(pi = pi, stats = stats))
  
  stats <- rep(0, 2)
  names(stats) <- c("lik", "df")
  
  for(gr in 1:ngroups){
    # gr = 1
    # cat(gr, fill = TRUE)
    
    fit_gr <- dm_fitOneGeneOneGroup(y = y[, igroups[[gr]], drop = FALSE], gamma0 = gamma0, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)
    
    if(any(is.na(fit_gr$pi))){
      pi <- matrix(NA, nrow = k, ncol = ngroups, dimnames = list(rownames(y), lgroups))
      stats <- rep(NA, 2)
      names(stats) <- c("lik", "df")
      return(list(pi = pi, stats = stats))  
    }
    
    pi[,gr] <- fit_gr$pi
    stats <- stats + fit_gr$stats
    
  }
  
  return(list(pi = pi, stats = stats))
  ### pi and stats can have NAs
}


