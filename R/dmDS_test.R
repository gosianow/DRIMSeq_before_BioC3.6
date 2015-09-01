#######################################################
#  group testing
#######################################################


# stats_full = x@fit_full@metadata; stats_null = x@fit_null@metadata

dmDS_test <- function(stats_full, stats_null){

  ## calculate lr
  cat("* Calculating likelihood ratio statistics.. \n")
  time_start <- Sys.time()
    
  lr <- 2*(rowSums(stats_full) - stats_null[, "lik"])
  
  nrgroups <- ncol(stats_full)
  
  df <- (nrgroups - 1)*stats_null[, "df"]
  
  pvalue <- pchisq(lr, df = df , lower.tail = FALSE)
  
  adj_pvalue <- p.adjust(pvalue, method="BH")
  
  table <- data.frame(gene_id = rownames(stats_full), lr = lr, df = df, pvalue = pvalue, adj_pvalue = adj_pvalue, stringsAsFactors = FALSE)

  o <- order(table[, "pvalue"])
  table <- table[o,]

  time_end <- Sys.time()
  cat("Took ", as.numeric(time_end - time_start), " seconds.\n")

  return(table)
  
}

