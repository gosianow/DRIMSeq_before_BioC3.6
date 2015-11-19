#' @include dm_cpm.R
NULL

# counts = x@counts; genotypes = x@genotypes; blocks = x@blocks; samples = x@samples; min_samps_gene_expr = 70; min_gene_expr = 1; min_samps_feature_prop = 5; min_feature_prop = 0.1; max_features = Inf; minor_allele_freq = 10; BPPARAM = BiocParallel::MulticoreParam(workers = 5)

dmSQTL_filter <- function(counts, genotypes, blocks, samples, min_samps_gene_expr = 70, min_gene_expr = 1, min_samps_feature_prop = 5, min_feature_prop = 0.1, max_features = Inf, minor_allele_freq = 5, BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  ########################################################
  # filtering on counts, put NA for samples with low gene expression
  ########################################################
  
  ### calculate cpm
  counts_for_cpm <- counts@unlistData
  counts_for_cpm[is.na(counts_for_cpm)] <- 0
  
  counts_cpm <- new("MatrixList", unlistData = dm_cpm(counts_for_cpm), partitioning = counts@partitioning) ### cpm can not handle NAs, so repalce NAs with 0s
  
  inds <- which(width(counts) > 1)

  counts_new <- lapply(inds, function(g){
    # g = 1
    expr_cpm_gene <- counts_cpm[[g]]
    expr_gene <- counts[[g]]
    
    ### genes with min expression
    if(! sum(colSums(expr_cpm_gene) >= min_gene_expr, na.rm = TRUE) >= min_samps_gene_expr )
      return(NULL)
      
      samps2keep <- colSums(expr_cpm_gene) > 0 & !is.na(expr_cpm_gene[1, ])
      
      if(sum(samps2keep) == 0)
        return(NULL)
    
    samps2keep <- colSums(expr_cpm_gene) >= min_gene_expr & !is.na(expr_cpm_gene[1, ])
    
    if(sum(samps2keep) < max(1, min_samps_feature_prop))
      return(NULL)
    
    prop <- prop.table(expr_gene[, samps2keep, drop = FALSE], 2) 
    trans2keep <- rowSums(prop >= min_feature_prop) >= min_samps_feature_prop
    
    ### no genes with one transcript
    if(sum(trans2keep) <= 1)
      return(NULL)
    
    #### Have to think how to order the transcripts because here I do not have the same grouping 
    # if(!max_features == Inf){
    #   if(sum(trans2keep) > max_features){
    #     tr_order <- order(-rowQuantiles(-prop, min_samps_feature_prop/ncol(prop)), decreasing = TRUE)
    #     trans2keep <- trans2keep[trans2keep]
    #     trans2keep <- names(trans2keep[sort(tr_order[1:max_features])])
    #   }
    # }
    
    expr <- expr_gene[trans2keep, , drop = FALSE] 
    expr[, !samps2keep] <- NA
    
    return(expr)
    
  })
  
  names(counts_new) <- names(counts)[inds]
  NULLs <- !sapply(counts_new, is.null)
  counts_new <- counts_new[NULLs]
  counts_new <- MatrixList(counts_new)
  
  ########################################################
  # filtering on genotypes
  ########################################################
  
  genotypes <- genotypes[inds[NULLs]]
  blocks <- blocks[inds[NULLs], ]
  
  genotypes_new <- BiocParallel::bplapply(1:length(counts_new), function(g){ 
    # g = 1
    
    counts_gene <- counts_new[[g]]
    genotypes_gene <- genotypes[[g]]
    
    ## NA for samples with non expressed genes and missing genotype
    genotypes_gene[, is.na(counts_gene[1,])] <- NA
    genotypes_gene[genotypes_gene == -1] <- NA
    
    ##### Keep genotypes with at least minor_allele_freq number of variants per group; in other case replace them with NAs
    genotypes_gene <- apply(genotypes_gene, 1, function(x){
      # x <- genotypes_gene[6,]
      
      tt <- table(x)
      
      if( length(tt)==1 )
        return(NULL)
      if( length(tt)==2 ){
        if(any(tt <= minor_allele_freq))
          return(NULL)
        return(x)
      }else{
        if(sum(tt <= minor_allele_freq) >= 2)
          return(NULL)
        x[x == names(tt[tt <= minor_allele_freq])] <- NA
        return(x)
      }    
    })
    
    if(!is.null(genotypes_gene)){
      if(is.list(genotypes_gene))
        genotypes_gene <- do.call(rbind, genotypes_gene)
      else
        genotypes_gene <- t(genotypes_gene)
    }
    
    return(genotypes_gene)
    
  }, BPPARAM = BPPARAM)
  
  names(genotypes_new) <- names(genotypes)
  NULLs <- !sapply(genotypes_new, is.null)
  genotypes_new <- genotypes_new[NULLs]
  genotypes_new <- MatrixList(genotypes_new)
  counts_new <- counts_new[NULLs]
  
  blocks <- blocks[NULLs, ]
  
  ########################################################
  # filtering on blocks
  ########################################################
  
  inds <- 1:length(genotypes_new)
  
  blocks_new <- MatrixList(lapply(inds, function(b){
    # b = 1
    blocks[[b]][blocks[[b]][, "block_id"] %in% rownames(genotypes_new[[b]]), , drop = FALSE]
  }))
  names(blocks_new) <- names(genotypes_new)

  data <- new("dmSQTLdata", counts = counts_new, genotypes = genotypes_new, blocks = blocks_new, samples = samples)
  
  return(data)
  
}









