#######################################################
#  group testing
#######################################################


# stats_full = x@fit_full@metadata; stats_null = x@fit_null[[1]]@metadata

dmDS_test <- function(stats_full, stats_null, test = "lr", n, verbose = FALSE){

  ## calculate lr
  if(verbose) cat("* Calculating likelihood ratio statistics.. \n")
  time_start <- Sys.time()
  
  lik_full <- stats_full[, grepl("lik_", colnames(stats_full)), drop = FALSE]
  lik_null <- stats_null[, "lik"]
  
  dev_full <- stats_full[, grepl("dev_", colnames(stats_full)), drop = FALSE]
  dev_null <- stats_null[, "dev"]
  
  nrgroups <- rowSums(!is.na(lik_full))
  
  switch(test, 
    
    lr = {
      
      lr <- 2 * (rowSums(lik_full) - lik_null)
      
      # lr <- dev_null - rowSums(dev_full) ## based on deviance
      
      df <- (nrgroups - 1) * stats_null[, "df"]
      
      df[nrgroups == 0] <- NA 
      lr[nrgroups == 0] <- NA 
      
      pvalue <- pchisq(lr, df = df , lower.tail = FALSE)
      
      adj_pvalue <- p.adjust(pvalue, method="BH")
      
      table <- data.frame(gene_id = rownames(stats_full), lr = lr, df = df, pvalue = pvalue, adj_pvalue = adj_pvalue, stringsAsFactors = FALSE)
      
      },
      
      f = {
        
        lr <- dev_null - rowSums(dev_full)
        
        df1 <- (nrgroups - 1) * stats_null[, "df"]
        df2 <- n - nrgroups * stats_null[, "df"]
        
        qdisp <- rowSums(dev_full) / df2
        
        lr[nrgroups == 0] <- NA 
        df1[nrgroups == 0] <- NA 
        qdisp[nrgroups == 0] <- NA
        
        
        f = (lr / df1) / qdisp
        
        pvalue <- pf(f, df1 = df1, df2 = df2, lower.tail = FALSE, log.p = FALSE)
        
        adj_pvalue <- p.adjust(pvalue, method="BH")
        
        table <- data.frame(gene_id = rownames(stats_full), f = f, df1 = df1, df2 = df2, pvalue = pvalue, adj_pvalue = adj_pvalue, stringsAsFactors = FALSE)
        
      }

    )
  

  # o <- order(table[, "pvalue"])
  # table <- table[o,]
  
  rownames(table) <- NULL

  time_end <- Sys.time()
  if(verbose) cat("Took ", as.numeric(time_end - time_start), " seconds.\n")

  return(table)
  
}

