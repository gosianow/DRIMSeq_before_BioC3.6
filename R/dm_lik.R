##############################################################################
# log-likelihood for k-1 parameters (pi)
##############################################################################

# pi <- pi[-length(pi)]

dm_lik <- function(pi, gamma0, y){
  ## pi has length of k-1
  ## gamma0 has length 1
  ## y has k rows and any number of columns
  
  k <- nrow(y)
  N <- ncol(y)
  S <- colSums(y)  
  l <- 0

  pi <- c(pi, 1 - sum(pi))
  
  for(j in 1:N){  
    # j=1
    l <- l - sum(log(gamma0 + 1:S[j] - 1))   
    for(i in 1:k){   
      # i=3
      if(y[i,j] == 0) lij <- 0
      else lij <- sum(log(pi[i] * gamma0 + 1:y[i,j] - 1))     
      l <- l + lij      
    }
  }
	
  return(l)
  
}

