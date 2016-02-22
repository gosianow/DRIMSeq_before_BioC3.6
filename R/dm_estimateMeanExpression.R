##############################################################################
# calculate common dispersion 
##############################################################################

dm_estimateMeanExpression <- function(counts, verbose = FALSE, 
  BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  ### calculate mean expression of genes 
  if(verbose) message("* Calculating mean gene expression.. \n")
  
  inds <- 1:length(counts)

  time <- system.time(mean_expression <- unlist(BiocParallel::bplapply(inds, function(g){ 

    mean(colSums(counts[[g]]), na.rm = TRUE)

    }, BPPARAM = BPPARAM)))
  
  names(mean_expression) <- names(counts)

  if(verbose) message("Took ", time["elapsed"], " seconds.\n")
  if(verbose) message("*** Mean gene expression: ", head(mean_expression), "... \n")
  
  
  return(mean_expression)

}


