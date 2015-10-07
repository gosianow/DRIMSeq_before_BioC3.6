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





##############################################################################
# Score-function for k-1 parameters (pi) 
##############################################################################

dm_score <- function(pi, gamma0, y){  

  k <- nrow(y)
  N <- ncol(y) 
  ykm1 <- y[-k, , drop=FALSE]
  yk <- y[k, ]

  pik <- 1-sum(pi)
  
  S <- rep(0, k-1)
  
  for(j in 1:N){ 
    # j=1
    if(yk[j] == 0) Skj <- 0
    else Skj <- sum(gamma0 / (pik * gamma0 + 1:yk[j] - 1))    
    
    for(i in 1:(k-1)){
      # i=1
      if(y[i,j] == 0) Sij <- 0
      else Sij <- sum(gamma0 / (pi[i] * gamma0 + 1:ykm1[i,j] - 1)) 
      S[i] <- S[i] + Sij
    }
    S <- S - Skj
  }
  
  return(S)
  
}

