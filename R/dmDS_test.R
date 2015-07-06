#######################################################
#  group testing
#######################################################


dmDS_test <- function(fit, verbose = FALSE, BPPARAM = MulticoreParam(workers=1)){
  
  ## calculate lr
  cat("Calculating likelihood ratio statistics.. \n")
  
  gene_list <- names(fit$fit_full)
	
  time <- system.time(lr_list <- bplapply(gene_list, function(g){
    # g = "FBgn0002121"
    
    if(verbose) cat("Testing gene: ", g, fill = TRUE)
    
    if(is.null(fit$fit_null[[g]]) || is.null(fit$fit_full[[g]])) 
    return(rep(NA, 6))
		
      lik_null <- fit$fit_null[[g]]$lik

      lik_full <- sum(fit$fit_full[[g]]$lik)

      lr <-  2*(lik_full - lik_null)
			
      nr_groups <- length(fit$fit_full[[g]]$df)
			
      # df <- DFfull - DFnull
      df <- fit$fit_null[[g]]$df * (nr_groups - 1) # (k-1) * nr of groups
			
      pvalue <- pchisq(lr, df = df , lower.tail = FALSE)
      
    return(c(lik_full, lik_null, nr_groups, df, lr, pvalue))
  }, BPPARAM = BPPARAM))
  
	cat("Took ", time["elapsed"], " seconds.\n")
	cat("Generating table with results.. \n")
	
  lr <- do.call(rbind, lr_list)
	colnames(lr) <- c("lik_full", "lik_null", "nr_groups", "df", "lr", "pvalue")
  adj_pvalue <- p.adjust(lr[, "pvalue"], method="BH")
  
	table <- data.frame(gene_id = gene_list, lr, adj_pvalue, stringsAsFactors = FALSE)
	
	o <- order(table[, "pvalue"])
	
  table <- table[o,]
  
  return(table)
  
  
}

