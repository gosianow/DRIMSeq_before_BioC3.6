##############################################################################
# log-likelihood for k parameters (pi)
##############################################################################
## Computes the log-likelihood -- with gamma functions -- 
dmLogLikG <- function(pi, gamma0, y){
  
  N <- ncol(y)
  S <- colSums(y)
  
  l <- sum( colSums( lgamma(y + repCol(pi, N) * gamma0) - lgamma(repCol(pi, N) * gamma0) ) )

  l <- N * lgamma(gamma0) - sum(lgamma(S + gamma0)) + l
  
  return(l)
  
}
