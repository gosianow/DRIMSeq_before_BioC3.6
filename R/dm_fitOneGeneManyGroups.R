### one gene, many groups


dm_fitOneGeneManyGroups <- function(y, ngroups, lgroups, igroups, gamma0, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE){
  # NULL for filtered genes or genes with one exon
  k <- dim(y)[1]
  if(k <= 1) return(NULL)
  
  pi <- matrix(0, nrow = k, ncol = ngroups)
	rownames(pi) <- rownames(y)
  lik = rep(NA, ngroups)
  df = rep(0, ngroups)
  
  for(gr in 1:ngroups){
    # gr = 1
    # cat(gr, fill = TRUE)
    
    fit_gr <- dm_fitOneGeneOneGroup(y = y[, igroups[[gr]], drop = FALSE], gamma0 = gamma0, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)
    
    if(is.null(fit_gr)) return(NULL)
    
    pi[,gr] <- fit_gr$pi
    lik[gr] <- fit_gr$lik
    df[gr] <- fit_gr$df
    
  }
  
  colnames(pi) <- names(df) <- names(lik) <- lgroups

	if(verbose){
		cat("gamma0", gamma0, "\n")
		cat("lik", lik, "\n")
	}
	
  return(list(pi = pi, gamma0 = gamma0, lik = lik, df = df))
  
}


