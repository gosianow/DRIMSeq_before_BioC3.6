##############################################################################
# log-likelihood for k-1 parameters (pi)
##############################################################################

# pi <- pi[-length(pi)]
## Computes the log-likelihood -- with gamma functions -- for k-1 parameters
dmLogLikGkm1 <- function(pi, gamma0, y){
  ## pi has length of k-1
  
  N <- ncol(y)
  S <- colSums(y)
  
  pi <- c(pi, 1 - sum(pi))
  # names(pi) <- rownames(y)
  
  l <- sum( colSums( lgamma(y + pi * gamma0) - lgamma(pi * gamma0) ) )
  
  l <- N * lgamma(gamma0) - sum(lgamma(S + gamma0)) + l
  
  return(l)
  
}

