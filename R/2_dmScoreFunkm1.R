##############################################################################
# Score-function for k-1 parameters (pi) 
##############################################################################

# pi <- pi[-length(pi)]

dmScoreFunkm1 <- function(pi, gamma0, y){  

  k <- nrow(y)
  N <- ncol(y) 
  ykm1 <- y[-k, , drop=FALSE]
  yk <- y[k, ] # as.numeric(y[k, ]) ## problem in 1:yk[j]

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

