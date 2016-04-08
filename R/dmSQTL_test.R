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




dmSQTL_permutations_all_genes <- function(x, fit_null, results, max_nr_perm_cycles = 10, max_nr_min_nr_sign_pval = 1e3, prop_mode, prop_tol, n, test, verbose, BPPARAM){
  
  fit_full <- x@fit_full
  
  nr_perm_tot <- 0
  nr_perm_cycles <- 0
  min_nr_sign_pval <- 0
  n <- ncol(x@counts)
  
  
  pval <- results$pvalue
  nas <- is.na(pval)
  pval <- pval[!nas]
  pval <- factor(pval)
  sum_sign_pval <- rep(0, length(pval))
  
  pval_perm_all <- matrix(NA, ncol = max_nr_perm_cycles, nrow = nrow(results))
  
  # ds_genes <- results$adj_pvalue < 0.1
  
  while(nr_perm_cycles < max_nr_perm_cycles && min_nr_sign_pval < max_nr_min_nr_sign_pval){

    permutation <- sample(n, n)
    
    ### Permute counts for all genes
    counts <- x@counts[, permutation, drop = FALSE]
    
    fit_full_perm <- dmSQTL_fitOneModel(counts = counts, genotypes = x@genotypes, dispersion = slot(x, x@dispersion), model = "full", prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
    
    fit_null_perm <- dmSQTL_fitOneModel(counts = counts, genotypes = x@genotypes, dispersion = slot(x, x@dispersion), model = "null", prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
    
    results_perm <- dmSQTL_test(fit_full = fit_full_perm, fit_null = fit_null_perm, test = test, n = n, verbose = verbose, BPPARAM = BPPARAM)
    
    nr_perm <- nrow(results_perm)
    nr_perm_tot <- nr_perm_tot + nr_perm
    nr_perm_cycles <- nr_perm_cycles + 1
    
    
    ### Count how many pval_permuted is lower than pval from the model
    pval_perm_all[, nr_perm_cycles] <- results_perm$pvalue
    nas_perm <- is.na(results_perm$pvalue)
    pval_perm <- results_perm$pvalue[!nas_perm]
    pval_perm_cut <- cut(pval_perm, c(-1, levels(pval), 2), right=FALSE)
    pval_perm_sum <- table(pval_perm_cut)
    pval_perm_cumsum <- cumsum(pval_perm_sum)[-length(pval_perm_sum)]
    names(pval_perm_cumsum) <- levels(pval)
    sum_sign_pval <- sum_sign_pval + pval_perm_cumsum[pval]
    
    pval_adj <- (sum_sign_pval + 1) / (nr_perm_tot + 1)
    
    min_nr_sign_pval <- min(sum_sign_pval)
    
    
    }
  
  
  pval_out <- rep(NA, nrow(results))
  pval_out[!nas] <- pval_adj
  
  return(list(pvalues_adjusted = pval_out, pvalues_permutations = pval_perm_all))
  
}












