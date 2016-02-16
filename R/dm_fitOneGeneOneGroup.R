##############################################################################
## estimate pi for given dispersion
##############################################################################


dm_fitOneGeneOneGroup <- function(y, gamma0, prop_mode = c("constrOptim", "constrOptimG")[2], prop_tol = 1e-12, verbose = FALSE){
  ### y must be features vs. samples
  ### If something is wrong, return NAs
  
  # NAs for genes with one feature
  kk <- nrow(y)
  if(kk < 2 || is.na(gamma0)) 
    return(list(pi = rep(NA, kk), stats = c(lik = NA, df = NA, dev = NA)))
  
  ### check for 0s in rows (features)
  keep_row <- rowSums(y) > 0
  if(sum(keep_row) < 2) 
    return(list(pi = rep(NA, kk), stats = c(lik = NA, df = NA, dev = NA))) ## must be at least two features
  
  y <- y[keep_row, , drop=FALSE]
  
  ### check for 0s in columns (replicates)
  keep_col <- colSums(y) > 0
  y <- y[, keep_col, drop=FALSE]
  
  pi_init <- rowSums(y)/sum(y)
  
  k <- length(pi_init) ## k - number of features
  
  if(sum(keep_col) == 1){
    
    keep_row[keep_row] <- pi_init
    pi <- keep_row
    
    df <- k - 1
    lik <- dm_likG(pi = pi_init[-k], gamma0 = gamma0, y = y) # likelihood
    dev <- dm_devG(pi = pi_init[-k], gamma0 = gamma0, y = y) # deviance
    
    stats <- c(lik, df, dev)
    names(stats) <- c("lik", "df", "dev")

    return(list(pi = pi, stats = stats))
  }

  ### Use more strict constrains (based on observed data) for the proportions
  # pi_init_m <- y/matrix(colSums(y), nrow = nrow(y), ncol = ncol(y), byrow = TRUE)
  # pi_init <- rowMeans(pi_init_m)
  # pi_init <- pi_init/sum(pi_init)

  # pi_init_minmax <- apply(pi_init_m, MARGIN = 1, function(x){  
  #   c(min(x, na.rm = TRUE), max(x, na.rm = TRUE))
  #   })
  
  # pi_init_minmax[1, ] <- pi_init_minmax[1, ] - 0.001
  # pi_init_minmax[1, pi_init_minmax[1, ] < 0] <- 0
  # pi_init_minmax[2, ] <- pi_init_minmax[2, ] + 0.001
  # pi_init_minmax[2, pi_init_minmax[2, ] > 1] <- 1
  
  # ### Make sure that the initial value is within the interval
  # # print(all(pi_init > pi_init_minmax[1, ] & pi_init < pi_init_minmax[2, ]))

  # if(!all(pi_init > pi_init_minmax[1, ])){
    
  #   indx_replace <- pi_init <= pi_init_minmax[1, ]
  #   # print(indx_replace)
  #   pi_init_minmax[1, indx_replace] <- pi_init[indx_replace] - 0.001
  #   pi_init_minmax[1, pi_init_minmax[1, ] < 0] <- 0
    
  # }
  
  # if(!all(pi_init < pi_init_minmax[2, ])){
    
  #   indx_replace <- pi_init >= pi_init_minmax[2, ]
  #   # print(indx_replace)
  #   pi_init_minmax[2, indx_replace] <- pi_init[indx_replace] + 0.001
  #   pi_init_minmax[2, pi_init_minmax[2, ] > 1] <- 1
    
  # }
  
  
  switch(prop_mode, 
         
         ### must have constraint for SUM pi = 1 --> sum(pi) < 1 + eps & sum(pi) > 1 - eps
         constrOptim = { ## for k-1 parameters
           # if(verbose) cat("\n gene:", colnames(y)[1], "gamma0:", gamma0, fill = TRUE)
           
           ui <- rbind(diag(rep(1, k-1), k-1), diag(rep(-1, k-1), k-1), rep(-1, k-1))
           ci <- c(rep(0, k-1), rep(-1, k-1), -1 + .Machine$double.eps) 

           co <- constrOptim(pi_init[-k], f = dm_lik, grad = dm_score, ui = ui, ci = ci, control = list(fnscale = -1, reltol = prop_tol), gamma0 = gamma0, y = y)
           
           pi <- co$par
           dev <- dm_devG(pi = pi, gamma0 = gamma0, y = y) 
           pi <- c(pi, 1-sum(pi))
           lik <- co$value
           
         }, 
         
         constrOptimG = { ## for k-1 parameters with Gamma functions
           # if(verbose) cat("\n gene:", colnames(y)[1], "gamma0:", gamma0, fill = TRUE)
           
           ui <- rbind(diag(rep(1, k-1), k-1), diag(rep(-1, k-1), k-1), rep(-1, k-1))
           ci <- c(rep(0, k-1), rep(-1, k-1), -1 + .Machine$double.eps) 
           
           # ui <- rbind(diag(rep(1, k-1), k-1), diag(rep(-1, k-1), k-1), rep(-1, k-1))
           # ci <- c(pi_init_minmax[1, 1:(k-1)], - pi_init_minmax[2, 1:(k-1)], -1 + .Machine$double.eps) 

           co <- constrOptim(pi_init[-k], f = dm_likG, grad = dm_scoreG, ui = ui, ci = ci, control = list(fnscale = -1, reltol = prop_tol), gamma0 = gamma0, y = y)
           
           pi <- co$par
           dev <- dm_devG(pi = pi, gamma0 = gamma0, y = y)
           pi <- c(pi, 1-sum(pi))
           lik <- co$value
           
         })
  
  keep_row[keep_row] <- pi
  pi <- keep_row
  
  df <- k - 1
  
  stats <- c(lik, df, dev)
  names(stats) <- c("lik", "df", "dev")
  
  return(list(pi = pi, stats = stats))
  
}







