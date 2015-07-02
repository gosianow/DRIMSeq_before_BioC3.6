##############################################################################
# calculate common dispersion 
##############################################################################

dm_estimateMeanExpression <- function(data, BPPARAM = MulticoreParam(workers=1)){
	
  ### calculate mean expression of genes 
	cat("Calculating mean gene expression.. \n")
	
  time <- system.time(mean_expression <- unlist(bplapply(data$counts, function(g){ mean(colSums(g), na.rm = TRUE) }, BPPARAM = BPPARAM)))
	
	cat("Took ", time["elapsed"], " seconds.\n")
	
	return(mean_expression)
	
	}