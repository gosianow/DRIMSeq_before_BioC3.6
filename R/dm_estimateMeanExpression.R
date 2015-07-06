##############################################################################
# calculate common dispersion 
##############################################################################

dm_estimateMeanExpression <- function(data, BPPARAM = MulticoreParam(workers = 1)){
	
  ### calculate mean expression of genes 
	cat("Calculating mean gene expression.. \n")
	
  time <- system.time(mean_expression <- do.call(rbind, (bplapply(data@counts, function(g){ 

  	c(mean(colSums(g), na.rm = TRUE), nrow(g))

  	}, BPPARAM = BPPARAM))))
  
  colnames(mean_expression) <- c("mean_expression", "nr_features")
	
	cat("Took ", time["elapsed"], " seconds.\n")
	
	return(mean_expression)
	
	}