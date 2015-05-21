##############################################################################
# Score-function for k-1 parameters (pi) 
##############################################################################

# pi <- pi[-length(pi)]

# dmScoreFunGkm1(pi, gamma0, y)

## Score-function -- with gamma functions -- for k-1 parameters
dmScoreFunGkm1 <- function(pi, gamma0, y){  
  ## pi has length of k-1
  
  k <- nrow(y)
  N <- ncol(y)
  ykm1 <- y[-k, , drop=FALSE]
  yk <- y[k,] 
  pik <- 1-sum(pi)
 
	S <- gamma0 * rowSums( digamma(ykm1 + pi * gamma0) - digamma(pi * gamma0) - matrix(digamma(yk + gamma0 * pik) - digamma(gamma0 * pik), nrow = k-1, ncol = N, byrow = TRUE) ) 
		
  ### versions with repCol 
	# yk <- as.numeric(y[k,]) ## was a problem if only y[k,]  
	# S <- gamma0 * rowSums( digamma(ykm1 + repCol(pi, N) * gamma0) - repCol(digamma(pi * gamma0), N) - repRow(digamma(yk + gamma0 * pik) - digamma(gamma0 * pik), k-1) ) 
	

  return(S)
  
} 



# repCol(digamma(pi * gamma0), N)

# matrix(digamma(pi * gamma0), nrow = k-1, ncol = N, byrow = FALSE)


# repRow(digamma(yk + gamma0 * pik) - digamma(gamma0 * pik), k-1)

# matrix(digamma(yk + gamma0 * pik) - digamma(gamma0 * pik), nrow = k-1, ncol = N, byrow = TRUE)



