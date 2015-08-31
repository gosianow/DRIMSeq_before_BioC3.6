##############################################################################
# calculate common dispersion 
##############################################################################

dm_estimateMeanExpression <- function(counts, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
	
  ### calculate mean expression of genes 
  cat("* Calculating mean gene expression.. \n")
  
  inds <- 1:length(counts)

  time <- system.time(mean_expression <- unlist(BiocParallel::bplapply(inds, function(g){ 

  	mean(colSums(counts[[g]]), na.rm = TRUE)

  	}, BPPARAM = BPPARAM)))
  
  names(mean_expression) <- names(counts)

  cat("Took ", time["elapsed"], " seconds.\n")
  if(verbose) cat("*** Mean gene expression: ", head(mean_expression), "... \n")
  
  
  return(mean_expression)

}


