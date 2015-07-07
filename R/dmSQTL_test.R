#######################################################
#  group testing
#######################################################

dmSQTL_test <- function(fit, verbose=FALSE, BPPARAM = MulticoreParam(workers=1)){
  
  ## calculate lr
  cat("Calculating likelihood ratio statistics.. \n")
  
  gene_list <- names(fit$fit_full)
  
  time <- system.time(lr_list <- bplapply(gene_list, function(g){
    # g = "ENSG00000188822.6"
    
    fit_full_g <- fit$fit_full[[g]]
    fit_null_g <- fit$fit_null[[g]]
    
    out_test <- matrix(NA, length(fit_full_g), 6)
    colnames(out_test) <- c("lik_full", "lik_null", "nr_groups", "df", "lr", "pvalue")
    
    for(i in 1:length(fit_full_g)){
      # i = 1
      if(is.null(fit_full_g[[i]]) || is.null(fit_null_g[[i]])) 
         next # out_test[i, ] <- NA
      
      out_test[i, "lik_full"] <- sum(fit_full_g[[i]]$logLik)
      out_test[i, "lik_null"] <- fit_null_g[[i]]$logLik
      out_test[i, "nr_groups"] <- length(fit_full_g[[i]]$df)
      out_test[i, "df"] <- fit_null_g[[i]]$df * ( out_test[i, "nr_groups"] - 1 )

    }
    
    NAs <- complete.cases(out_test[, c("lik_full", "lik_null", "nr_groups", "df"), drop = FALSE])
    
    # lr <-  2 * (lik_full - lik_null)
    out_test[NAs, "lr"] <- 2 * (out_test[NAs, "lik_full", drop = FALSE] - out_test[NAs, "lik_null", drop = FALSE])
    out_test[NAs, "pvalue"] <- pchisq(out_test[NAs, "lr", drop = FALSE], df = out_test[NAs, "df", drop = FALSE] , lower.tail = FALSE)
    
    return(out_test)
    
  }, BPPARAM = BPPARAM))
	
	cat("Took ", time["elapsed"], " seconds.\n")
	cat("Generating table with results.. \n")
	
	### gene and snp IDs
  id_list <- bplapply(gene_list, function(g){
    # g = gene_list[1] 
		
		matrix(c(rep(g, length(fit$fit_full[[g]])), names(fit$fit_full[[g]])), nrow = length(fit$fit_full[[g]]), byrow = FALSE, dimnames = list(NULL, c("gene_id", "snp_id")))
	 
  }, BPPARAM = BPPARAM)
  

	id <- do.call(rbind, id_list)	
  lr <- data.matrix(do.call(rbind, lr_list))
  adj_pvalue <- p.adjust(lr[, "pvalue"], method="BH")
	
	table <- data.frame(id, lr, adj_pvalue, stringsAsFactors = FALSE)
	
  o <- order(table[, "pvalue"])  
  
	table <- table[o,]


  return(table)
  
  
}


