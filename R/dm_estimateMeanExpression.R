##############################################################################
# calculate common dispersion 
##############################################################################

dm_estimateMeanExpression <- function(counts, BPPARAM = MulticoreParam(workers = 1)){
	
  ### calculate mean expression of genes 
  cat("Calculating mean gene expression.. \n")

  time <- system.time(mean_expression <- unlist(BiocParallel::bplapply(counts, function(g){ 

  	mean(colSums(g), na.rm = TRUE)

  	}, BPPARAM = BPPARAM)))


  cat("Took ", time["elapsed"], " seconds.\n")
  cat("** Mean gene expression: ", head(mean_expression), "... \n")
  
  
  return(mean_expression)

}


