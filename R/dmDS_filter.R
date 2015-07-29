# data must be dmDSdata object
# min_samps_gene_expr = 3; min_gene_expr = 1; min_samps_feature_prop = 3; min_feature_prop = 0.01; max_features = Inf


dmDS_filter <- function(counts, samples, min_samps_gene_expr = 3, min_gene_expr = 1, min_samps_feature_prop = 3, min_feature_prop = 0.01, max_features = Inf){
  
  ### calculate cpm
  expr_cpm_spl <- new("MatrixList", unlistData = edgeR::cpm(counts@unlistData), partitioning = counts@partitioning)
  
  expr_spl <- counts 
  
  counts <- lapply(names(expr_spl), function(g){
    # g = "FBgn0000008"
    # print(g)
    expr_cpm_gene <- expr_cpm_spl[[g]]
    expr_gene <- expr_spl[[g]]
    
    ### no genes with one transcript
    if(dim(expr_gene)[1] == 1)
      return(NULL)
    
    ### genes with min expression
    if(! sum(colSums(expr_cpm_gene) > min_gene_expr) >= min_samps_gene_expr )
      return(NULL)
    
    samps2keep <- colSums(expr_cpm_gene) != 0 & !is.na(expr_cpm_gene[1, ])
    
    prop <- prop.table(expr_gene[, samps2keep], 2) 
    trans2keep <- rowSums(prop > min_feature_prop) >= min_samps_feature_prop
    
    ### no genes with one transcript
    if(sum(trans2keep) <= 1)
      return(NULL)
    
    
    if(!max_features == Inf){
      if(sum(trans2keep) > max_features){
        
        tr_order <- order(apply(aggregate(t(prop[trans2keep, ]), by = list(group = samples$group[samps2keep]), median)[, -1], 2, max), decreasing = TRUE)
        
        trans2keep <- trans2keep[trans2keep]
        
        trans2keep <- names(trans2keep[sort(tr_order[1:max_features])])
        
      }
    }
    
    expr <- expr_gene[trans2keep, ] 
    
    return(expr)
    
  })
  
  names(counts) <- names(expr_spl)

  data_filtered <- new("dmDSdata", counts = MatrixList(counts), samples = samples)

  return(data_filtered)
  
}













