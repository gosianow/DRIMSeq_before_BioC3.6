#######################################################
#  group testing
#######################################################
# fit_full = dt@fit_full; fit_null = dt@fit_null; BPPARAM = BiocParallel::MulticoreParam(workers = 10)


dmSQTL_test <- function(fit_full, fit_null, test = "lr", n, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  ## calculate lr
  if(verbose) cat("* Calculating likelihood ratio statistics.. \n")
  time_start <- Sys.time()
  
  inds <- 1:length(fit_full)
  gene_list <- names(fit_full)
  
  switch(test, 
    
    lr = {
      
      table_list <- BiocParallel::bplapply(inds, function(g){
        # g = 662
        
        lik_full <- fit_full[[g]]@metadata[, grepl("lik_", colnames(fit_full[[g]]@metadata)), drop = FALSE]
        lik_null <- fit_null[[g]]@metadata[, "lik"]
        
        lr <- 2*(rowSums(lik_full, na.rm = TRUE) - lik_null)
        
        nrgroups <- rowSums(!is.na(lik_full))
        
        df <- (nrgroups - 1)*fit_null[[g]]@metadata[, "df"] ### negative when NAs in all groups in lik
        
        df[nrgroups == 0] <- NA 
        lr[nrgroups == 0] <- NA 
        
        pvalue <- pchisq(lr, df = df , lower.tail = FALSE)
        
        tt <- data.frame(gene_id = gene_list[g], snp_id = rownames(fit_full[[g]]@metadata), lr = lr, df = df, pvalue = pvalue, stringsAsFactors = FALSE)
        
        return(tt)
        
        }, BPPARAM = BPPARAM)
      
    },
    
    fql = {

      table_list <- BiocParallel::bplapply(inds, function(g){
        # g = 662
        
        dev_full <- fit_full[[g]]@metadata[, grepl("dev_", colnames(fit_full[[g]]@metadata)), drop = FALSE]
        dev_null <- fit_null[[g]]@metadata[, "dev"]

        # lr <- dev_null - rowSums(dev_full, na.rm = TRUE) ### This gives negative values!

        lik_full <- fit_full[[g]]@metadata[, grepl("lik_", colnames(fit_full[[g]]@metadata)), drop = FALSE]
        lik_null <- fit_null[[g]]@metadata[, "lik"]
        
        lr <- 2*(rowSums(lik_full, na.rm = TRUE) - lik_null)
      
        nrgroups <- rowSums(!is.na(lik_full))
        
        df1 <- (nrgroups - 1) * fit_null[[g]]@metadata[, "df"]
        df2 <- (n[g] - nrgroups) * fit_null[[g]]@metadata[, "df"] ### (n-J)*(q-1)
        
        qdisp <- rowSums(dev_full, na.rm = TRUE) 
        
        f = (lr / df1) / (qdisp / df2)
        
        pvalue <- pf(f, df1 = df1, df2 = df2, lower.tail = FALSE, log.p = FALSE)
        
        tt <- data.frame(gene_id = gene_list[g], snp_id = rownames(fit_full[[g]]@metadata), f = f, df1 = df1, df2 = df2, pvalue = pvalue, stringsAsFactors = FALSE)
        
        return(tt)
        
        }, BPPARAM = BPPARAM)
      
    },
    
    f = {

      table_list <- BiocParallel::bplapply(inds, function(g){
        # g = 662
        
        lik_full <- fit_full[[g]]@metadata[, grepl("lik_", colnames(fit_full[[g]]@metadata)), drop = FALSE]
        lik_null <- fit_null[[g]]@metadata[, "lik"]
        
        lr <- 2*(rowSums(lik_full, na.rm = TRUE) - lik_null)
      
        nrgroups <- rowSums(!is.na(lik_full))
        
        df1 <- (nrgroups - 1) * fit_null[[g]]@metadata[, "df"]
        df2 <- (n[g] - nrgroups) * fit_null[[g]]@metadata[, "df"] ### (n-J)*(q-1)
        
        f <- lr / df1
        
        pvalue <- pf(f, df1 = df1, df2 = df2, lower.tail = FALSE, log.p = FALSE)
        
        tt <- data.frame(gene_id = gene_list[g], snp_id = rownames(fit_full[[g]]@metadata), f = f, df1 = df1, df2 = df2, pvalue = pvalue, stringsAsFactors = FALSE)
        
        return(tt)
        
        }, BPPARAM = BPPARAM)
      
    }
    
    )
  

  table <- do.call(rbind, table_list)
  
  adj_pvalue <- p.adjust(table[, "pvalue"], method="BH")
  
  table$adj_pvalue <- adj_pvalue
  
  # o <- order(table[, "pvalue"])  
  # table <- table[o,]
  
  rownames(table) <- NULL
  
  time_end <- Sys.time()
  if(verbose) cat("Took ", as.numeric(time_end - time_start), " seconds.\n")
  
  return(table)
  
  
}


