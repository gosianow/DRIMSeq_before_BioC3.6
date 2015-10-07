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

##############################################################################
# observed FIM for k-1 parameters (pi) // FIXED mistakes
##############################################################################

# pi <- pi[-length(pi)]

#### Computing the inverce of observed Fisher Information Matrix
dm_FIMobs <- function(pi, gamma0, y, inv = TRUE){  
  ## pi has length of k-1

  k <- nrow(y)
  N <- ncol(y)
  D <- rep(0, k-1)
  Dil <- 0
  pik <- 1-sum(pi)
  
  for(j in 1:N){
    # j=1
    if(y[k, j] == 0) Dil <- Dil + 0
    else Dil <- Dil +  sum(-gamma0^2 / (pik * gamma0 + 1:y[k, j] - 1) ^2) 

    for(i in 1:(k-1)){
      # i=1
      if(y[i,j] == 0) Dii <- 0
      else Dii <- sum(-gamma0^2 / (pi[i] * gamma0 + 1:y[i,j] - 1) ^2) 
      D[i] <- D[i] + Dii
    }
  }
  
  H <- matrix(Dil, k-1, k-1)
  diag(H) <- diag(H) + D
  
  if(!inv)
    return(-H)
  
  invFIM <- solve(-H)
  return(invFIM)
}

##############################################################################
# expected FIM for k-1 parameters (pi) 
##############################################################################


#### Beta-Binomial
## does not work for big gamma0 !!!
dBetaBin <- function(x, n, a, b){
  beta(x+a, n-x+b) / beta(a,b) * choose(n,x)
}

pBetaBin <- function(x, n, a, b) {
  ## x must be >=1
  sapply(x, function(xx){ 
    1 - sum(dBetaBin(0:(xx-1), n, a, b))
    })
}

## works for big gamma0
dBetaBin2 <- function(x,n,a,b){ 
  sapply(x, function(xx){
    res <- lchoose(n,xx)
    if(xx==0) res <- res else res <- res + sum(log(a + 1:xx - 1))
    if(xx==n) res <- res else res <- res + sum(log(b + 1:(n-xx) -1))
    res <- res - sum(log(a + b + 1:n - 1))
    names(res) <- NULL
    exp(res)
  })
}

pBetaBin2 <- function(x, n, a, b) {
  ## x must be >=1
  sapply(x, function(xx){ 
    1 - sum(dBetaBin2(0:(xx-1), n, a, b))
  })
}


# x <- 1:5
# n <- 10
# a <- b <- 5
# 
# x = 1:y[i,k]; n = S[i]; a = pik * gamma0; b = (1 - pik) * gamma0

# dBetaBin(x, n, a, b)
# pBetaBin(x, n, a, b)
# 
# dBetaBin2(x,n,a,b)
# pBetaBin2(x,n,a,b)



  
#### Computing the inverce of expected Fisher Information Matrix
dm_FIMexp <- function(pi, gamma0, y, inv = TRUE){    
  ## pi has length of k-1
  
  k <- nrow(y)
  N <- ncol(y)
  D <- rep(0, k-1)
  Dil <- 0
  pik <- 1-sum(pi)
  S <- colSums(y)
  
  for(j in 1:N){
    # j=1
    if(y[k, j] == 0) Dil <- Dil + 0
    else Dil <- Dil +  sum(-gamma0^2 * pBetaBin2(1:y[k,j], n = S[j], a = pik * gamma0, b = (1 - pik) * gamma0) / (pik * gamma0 + 1:y[k, j] - 1) ^2) 
    
    for(i in 1:(k-1)){
      # i=1
      if(y[i,j] == 0) Dii <- 0
      else Dii <- sum(-gamma0^2 * pBetaBin2(1:y[i,j], n = S[j], a = pi[i]*gamma0, b = (1 - pi[i]) * gamma0) / (pi[i] * gamma0 + 1:y[i,j] - 1) ^2) 
      D[i] <- D[i] + Dii
    }
  }
  
  H <- matrix(Dil, k-1, k-1)
  diag(H) <- diag(H) + D
  
  if(!inv)
    return(-H)
  
  invFIM <- solve(-H)
  return(invFIM)
  
}




