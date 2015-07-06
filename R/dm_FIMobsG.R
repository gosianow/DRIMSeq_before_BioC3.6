##############################################################################
# observed FIM for k-1 parameters (pi) // FIXED mistakes
##############################################################################

# pi <- pi[-length(pi)]

### FIM -- with gamma functions -- for k-1 parameters
dm_FIMobsG <- function(pi, gamma0, y, inv = TRUE){  
  ## pi has length of k-1
  
  k <- nrow(y)
  N <- ncol(y)
  D <- rep(0, k-1)
  Dil <- 0
  pik <- 1-sum(pi)
  ykm1 <- y[-k, , drop=FALSE]
  
Dil <- gamma0^2 * sum(trigamma(y[k,] + gamma0 * pik) - trigamma(gamma0 * pik))
  
D <- gamma0^2 * rowSums(trigamma(ykm1 + pi * gamma0) - trigamma(pi * gamma0))

  H <- matrix(Dil, k-1, k-1)
  diag(H) <- diag(H) + D
  
  if(!inv)
    return(-H)
  
  invFIM <- solve(-H)
  return(invFIM)
}

