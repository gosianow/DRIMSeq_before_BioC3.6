##############################################################################
## Computes the log-likelihood -- with gamma functions -- for k-1 parameters
##############################################################################

# pi <- pi[-length(pi)]

dm_likG <- function(pi, gamma0, y){
  ## pi has length of k-1
  ## gamma0 has length 1
  ## y has k rows and any number of columns
  
  N <- ncol(y)
  S <- colSums(y)
  
  pi <- c(pi, 1 - sum(pi))
  
  l <- sum( colSums( lgamma(y + pi * gamma0) - lgamma(pi * gamma0) ) )
  
  l <- N * lgamma(gamma0) - sum(lgamma(S + gamma0)) + l
  
  return(l)
  
}

