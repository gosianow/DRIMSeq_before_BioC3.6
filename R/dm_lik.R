##############################################################################
# log-likelihood for k-1 parameters (pi)
##############################################################################

# pi <- pi[-length(pi)]

## Computes the log-likelihood
dm_lik <- function(pi, gamma0, y){
  ## pi has length of k-1
  
  k <- nrow(y)
  N <- ncol(y)
  S <- colSums(y)  
  l <- 0

  pi <- c(pi, 1 - sum(pi))
  names(pi) <- rownames(y)
  
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
	
	# cat("pi:", pi, fill = T)
  # cat("fn:", l, fill = T)
  return(l)
  
}

