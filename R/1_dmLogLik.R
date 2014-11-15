##############################################################################
# log-likelihood for k parameters (pi)
##############################################################################


## Computes the log-likelihood
dmLogLik <- function(pi, gamma0, y){
  
  k <- nrow(y) ## exons
  N <- ncol(y) ## samples
  S <- colSums(y) ## tot reads in samples
  l <- 0
  
  for(j in 1:N){  
    # j=1
    l <- l - sum(sapply(1:S[j], function(r){log(gamma0 + r - 1)}))    
    for(i in 1:k){   
      # i=3
      if(y[i,j] == 0) lij <- 0
      else lij <- sum(sapply(1:y[i,j], function(r){log(pi[i] * gamma0 + r - 1)}))    
      l <- l + lij      
    }
  }
  
  return(l)
  
}
