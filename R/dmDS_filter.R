
# counts = x@counts; samples = x@samples; min_samps_gene_expr = 6; min_gene_expr = 10; min_samps_feature_expr = 3; min_feature_expr = 10; min_samps_feature_prop = 3; min_feature_prop = 0.01; max_features = Inf

dmDS_filter <- function(counts, samples, min_samps_gene_expr = 6, min_gene_expr = 10, min_samps_feature_expr = 3, min_feature_expr = 10, min_samps_feature_prop = 3, min_feature_prop = 0.01, max_features = Inf){
  
  inds <- which(width(counts) > 1)
  
  counts_new <- lapply(inds, function(g){
    # g = 117
    # print(g)
    
    expr_features <- counts[[g]]
    
    ### genes with min expression
    if(! sum(colSums(expr_features) >= min_gene_expr, na.rm = TRUE) >= min_samps_gene_expr )
      return(NULL)
    
    ### features with min expression
    features2keep <- rowSums(expr_features >= min_feature_expr, na.rm = TRUE) >= min_samps_feature_expr
    
    ### no genes with one feature
    if(sum(features2keep) <= 1)
      return(NULL)
    
    expr_features <- expr_features[features2keep, , drop = FALSE]
    
    
    ### genes with zero expression
    samps2keep <- colSums(expr_features) > 0 & !is.na(expr_features[1, ])
    
    if(sum(samps2keep) < max(1, min_samps_feature_prop))
      return(NULL)
    
    
    prop <- prop.table(expr_features[, samps2keep, drop = FALSE], 2) 
    # prop.table(matrix(c(1,0), 2, 1), 2)
    # prop.table(matrix(c(0,0), 2, 1), 2)
    # prop.table(matrix(c(0,0, 1, 0), 2, 2), 2)
    
    ### features with min proportion
    features2keep <- rowSums(prop >= min_feature_prop) >= min_samps_feature_prop
    
    ### no genes with one feature
    if(sum(features2keep) <= 1)
      return(NULL)
    
    
    if(!max_features == Inf){
      if(sum(features2keep) > max_features){
        
        tr_order <- order(apply(aggregate(t(prop[features2keep, , drop = FALSE]), by = list(group = samples$group[samps2keep]), median)[, -1], 2, max), decreasing = TRUE)
        
        features2keep <- features2keep[features2keep]
        
        features2keep <- names(features2keep[sort(tr_order[1:max_features])])
        
      }
    }
    
    expr <- expr_features[features2keep, , drop = FALSE] 
    
    return(expr)
    
  })
  
  names(counts_new) <- names(counts)[inds]
  counts_new <- counts_new[!sapply(counts_new, is.null)]
  counts_new <- MatrixList(counts_new)
  
  data <- new("dmDSdata", counts = counts_new, samples = samples)
  
  return(data)
  
}













