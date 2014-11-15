
##############################################################################
# adjustements to profile likelihood
##############################################################################

# g=1
# piH <- fit[[g]]$piH
# y <- y[[g]]


dmAdjCROneGeneGroup <- function(y, gamma0, piH){
  # NULL for filtered genes or genes with one exon
  if(dim(y)[1] <= 1) return(NULL)
  
  ### y must be exons vs. samples
  
  ### check for 0s in rows (exons)
  keepRow <- rowSums(y) > 0
  if(sum(keepRow) < 2) return(NULL) ## must be at least two exons or transcripts
  y <- y[keepRow, , drop=FALSE]
  
  ### check for 0s in cols (samples)
  keepCol <- colSums(y) > 0
  if(sum(keepCol) < 2) return(NULL) ## must be at least two samples in a condition
  y <- y[, keepCol, drop=FALSE]
  
  N <- ncol(y) 
  piH <- piH[keepRow]
  
  ## 1/2 * log(det(N* FIM))
#   adjCR <- log(det(N * dmInvObsFIMkm1(pi = piH[-length(piH)], gamma0, y, inv = FALSE)))/2 ## if piH is NULL then returns -Inf
# if(adjCR == -Inf)
#   return(NULL)

  adjCR <- log(det(N * dmInvObsFIMGkm1(pi = piH[-length(piH)], gamma0, y, inv = FALSE) ))/2 
## with Gamma functions ## if piH is NULL then:
# Error in is.data.frame(x) :
#   dims [product 6] do not match the length of object [0]
  
  return(adjCR)
  
}


