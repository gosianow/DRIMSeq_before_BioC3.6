#######################################################
#  group testing
#######################################################

dmSQTL_test <- function(fit_full, fit_null, BPPARAM = MulticoreParam(workers=1)){
  
  ## calculate lr
  cat("Calculating likelihood ratio statistics.. \n")
  time_start <- Sys.time()
  
  gene_list <- names(fit_full)
  
  table_list <- BiocParallel::bplapply(gene_list, function(g){
    # g = "ENSG00000131037.8"
    
    lr <- 2*(fit_full[[g]]@statistics[, "lik"] - fit_null[[g]]@statistics[, "lik"])
    
    df <- fit_full[[g]]@statistics[, "df"] - fit_null[[g]]@statistics[, "df"]
    
    pvalue <- pchisq(lr, df = df , lower.tail = FALSE)
    
    tt <- S4Vectors::DataFrame(gene_id = g, snp_id = rownames(fit_full[[g]]@statistics), lr = lr, df = df, pvalue = pvalue, row.names = paste0(g, ":", rownames(fit_full[[g]]@statistics)))
    
    }, BPPARAM = BPPARAM)
  
  
  table <- do.call(rbind, table_list)
  
  adj_pvalue <- p.adjust(table[, "pvalue"], method="BH")
  
  table$adj_pvalue <- adj_pvalue
  
  o <- order(table[, "pvalue"])  
  table <- table[o,]
  
  time_end <- Sys.time()
  cat("Took ", as.numeric(time_end - time_start), " seconds.\n")
  
  return(table)
  
  
}


