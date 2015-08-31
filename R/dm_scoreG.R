##############################################################################
# Score-function for k-1 parameters (pi) 
##############################################################################

## Score-function -- with gamma functions -- for k-1 parameters
dm_scoreG <- function(pi, gamma0, y){  
  ## pi has length of k-1
  
  k <- nrow(y)
  N <- ncol(y)
  ykm1 <- y[-k, , drop=FALSE]
  yk <- y[k,] 
  pik <- 1-sum(pi)
 
	S <- gamma0 * rowSums( digamma(ykm1 + pi * gamma0) - digamma(pi * gamma0) - matrix(digamma(yk + gamma0 * pik) - digamma(gamma0 * pik), nrow = k-1, ncol = N, byrow = TRUE) ) 
		
  return(S)
  
} 



