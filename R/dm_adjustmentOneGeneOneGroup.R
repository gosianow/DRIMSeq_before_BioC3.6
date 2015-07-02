##############################################################################
# Cox-Reid adjustment to profile likelihood
##############################################################################


dm_adjustmentOneGeneOneGroup <- function(y, gamma0, pi){
  # NULL for filtered genes or genes with one exon
  if(dim(y)[1] <= 1) return(NULL)
  
  ### y must be exons vs. samples
  
  ### check for 0s in rows (exons)
  keep_row <- rowSums(y) > 0
  if(sum(keep_row) < 2) return(NULL) ## must be at least two exons or transcripts
  y <- y[keep_row, , drop=FALSE]
  
  ### check for 0s in cols (samples)
  keep_col <- colSums(y) > 0
  if(sum(keep_col) < 2) return(NULL) ## must be at least two samples in a condition
  y <- y[, keep_col, drop=FALSE]
  
  N <- ncol(y) 
  pi <- pi[keep_row]
  
  adj <- log(det(N * dm_FIMobsG(pi = pi[-length(pi)], gamma0, y, inv = FALSE) ))/2 
## with Gamma functions ## if pi is NULL then:
# Error in is.data.frame(x) :
#   dims [product 6] do not match the length of object [0]
  
  return(adj)
  
}


