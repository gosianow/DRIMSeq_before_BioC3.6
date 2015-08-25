#######################################################
#  group testing
#######################################################

dmSQTL_test <- function(fit_full, fit_null, BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  ## calculate lr
  cat("Calculating likelihood ratio statistics.. \n")
  time_start <- Sys.time()
  
  inds <- 1:length(fit_full)
  gene_list <- names(fit_full)
  
  table_list <- BiocParallel::bplapply(inds, function(g){
    # g = "ENSG00000131037.8"
    
    lr <- 2*(rowSums(fit_full[[g]]@metadata, na.rm = TRUE) - fit_null[[g]]@metadata[, "lik"])
    
    nrgroups <- apply(fit_full[[g]]@metadata, 1, function(r){ sum(!is.na(r)) })
    
    df <- (nrgroups - 1)*fit_null[[g]]@metadata[, "df"]
    
    pvalue <- pchisq(lr, df = df , lower.tail = FALSE)
    
    tt <- data.frame(gene_id = gene_list[g], snp_id = rownames(fit_full[[g]]@metadata), lr = lr, df = df, pvalue = pvalue, row.names = paste0(gene_list[g], ":", rownames(fit_full[[g]]@metadata)), stringsAsFactors = FALSE)
    
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


