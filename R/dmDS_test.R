#######################################################
#  group testing
#######################################################


# stats_full = x@statistics_full; stats_null = x@statistics_null

dmDS_test <- function(stats_full, stats_null){

  ## calculate lr
  cat("Calculating likelihood ratio statistics.. \n")
  time_start <- Sys.time()
    
  lr <- 2*(stats_full[, "lik"] - stats_null[, "lik"])
  
  df <- stats_full[, "df"] - stats_null[, "df"]
  
  pvalue <- pchisq(lr, df = df , lower.tail = FALSE)
  
  adj_pvalue <- p.adjust(pvalue, method="BH")
  
  table <- IRanges::DataFrame(gene_id = rownames(stats_full), lr = lr, df = df, pvalue = pvalue, adj_pvalue = adj_pvalue, row.names = rownames(stats_full))

  o <- order(table[, "pvalue"])
  table <- table[o,]

  time_end <- Sys.time()
  cat("Took ", as.numeric(time_end - time_start), " seconds.\n")

  return(table)
  
}

