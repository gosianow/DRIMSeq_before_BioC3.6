#' @include dm_cpm.R
NULL

# counts = x@counts; samples = x@samples; min_samps_gene_expr = 1; min_gene_expr = 1; min_samps_feature_prop = 1; min_feature_prop = 0.01; max_features = Inf

dmDS_filter <- function(counts, samples, min_samps_gene_expr = 3, min_gene_expr = 1, min_samps_feature_prop = 3, min_feature_prop = 0.01, max_features = Inf){
  
  ### calculate cpm
  counts_cpm <- new("MatrixList", unlistData = dm_cpm(counts@unlistData), partitioning = counts@partitioning)
  
  inds <- which(width(counts) > 1)
  
  counts_new <- lapply(inds, function(g){
    # g = 117
    # print(g)
    expr_cpm_gene <- counts_cpm[[g]]
    expr_gene <- counts[[g]]
    
    # ### no genes with one transcript
    # if(dim(expr_gene)[1] == 1)
    #   return(NULL)
    
    ### genes with min expression
    if(! sum(colSums(expr_cpm_gene) >= min_gene_expr) >= min_samps_gene_expr )
      return(NULL)
    
    samps2keep <- colSums(expr_cpm_gene) > 0 & !is.na(expr_cpm_gene[1, ])
    
    if(sum(samps2keep) == 0)
      return(NULL)
    
    prop <- prop.table(expr_gene[, samps2keep, drop = FALSE], 2) 
    # prop.table(matrix(c(1,0), 2, 1), 2)
    # prop.table(matrix(c(0,0), 2, 1), 2)
    # prop.table(matrix(c(0,0, 1, 0), 2, 2), 2)
    
    ### transcripts with min proportion
    trans2keep <- rowSums(prop >= min_feature_prop) >= min_samps_feature_prop
    
    ### no genes with one transcript
    if(sum(trans2keep) <= 1)
      return(NULL)
    
    
    if(!max_features == Inf){
      if(sum(trans2keep) > max_features){
        
        tr_order <- order(apply(aggregate(t(prop[trans2keep, , drop = FALSE]), by = list(group = samples$group[samps2keep]), median)[, -1], 2, max), decreasing = TRUE)
        
        trans2keep <- trans2keep[trans2keep]
        
        trans2keep <- names(trans2keep[sort(tr_order[1:max_features])])
        
      }
    }
    
    expr <- expr_gene[trans2keep, , drop = FALSE] 
    
    return(expr)
    
  })
  
  names(counts_new) <- names(counts)[inds]
  counts_new <- counts_new[!sapply(counts_new, is.null)]
  counts_new <- MatrixList(counts_new)
  
  data <- new("dmDSdata", counts = counts_new, samples = samples)

  return(data)
  
}













