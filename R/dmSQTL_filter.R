# counts = x@counts; genotypes = x@genotypes; samples = x@samples; min_samps_gene_expr = 70; min_gene_expr = 1; min_samps_feature_prop = 5; min_feature_prop = 0.1; max_features = Inf; minor_allel_freq = 0.10; BPPARAM = MulticoreParam(workers = 10)


dmSQTL_filter <- function(counts, genotypes, samples, min_samps_gene_expr = 70, min_gene_expr = 1, min_samps_feature_prop = 5, min_feature_prop = 0.1, max_features = Inf, minor_allel_freq = 0.05, BPPARAM = MulticoreParam(workers = 1)){

########################################################
# filtering on counts 
# put NA for samples with low gene expression
########################################################

  ### calculate cpm
  counts_for_cpm <- counts@unlistData
  counts_for_cpm[is.na(counts_for_cpm)] <- 0
  
  counts_cpm <- new("MatrixList", unlistData = edgeR::cpm(counts_for_cpm), partitioning = counts@partitioning) ### cpm can not handle NAs, so repalce NAs with 0s
  
  
  gene_list <- names(counts)
  
  counts_new <- lapply(gene_list, function(g){
    # g = "ENSG00000008130.10"
    # print(g)
    expr_cpm_gene <- counts_cpm[[g]]
    expr_gene <- counts[[g]]
    
    ### no genes with one transcript
    if(dim(expr_gene)[1] == 1)
    return(NULL)
    
    ### genes with min expression
    if(! sum(colSums(expr_cpm_gene) > min_gene_expr, na.rm = TRUE) >= min_samps_gene_expr )
    return(NULL)
    
    samps2keep <- colSums(expr_cpm_gene) > min_gene_expr & !is.na(expr_cpm_gene[1, ])
    
    prop <- prop.table(expr_gene[, samps2keep], 2) 
    trans2keep <- rowSums(prop > min_feature_prop) >= min_samps_feature_prop
    
    ### no genes with one transcript
    if(sum(trans2keep) <= 1)
    return(NULL)
    
    #### Have to think how to order the transcripts because here I do not have the same grouping 
    # if(!max_features == Inf){
    #   if(sum(trans2keep) > max_features){

    #     require("matrixStats")

    #     tr_order <- order(-rowQuantiles(-prop, min_samps_feature_prop/ncol(prop)), decreasing = TRUE)

    #     trans2keep <- trans2keep[trans2keep]

    #     trans2keep <- names(trans2keep[sort(tr_order[1:max_features])])

    #   }
    # }
    
    expr <- expr_gene[trans2keep, ] 
    expr[, !samps2keep] <- NA

    return(expr)
    
    })

names(counts_new) <- names(counts)
counts_new <- MatrixList(counts_new)

########################################################
# filtering on genotypes
########################################################

minor_allel_nr <- ceiling(minor_allel_freq * nrow(samples))

gene_list <- names(counts_new)

genotypes_new <- bplapply(gene_list, function(g){ 
  # g <- gene_list[1]; print(g)

  counts_gene <- counts[[g]]
  genotypes_gene <- genotypes[[g]]

  ## NA for samples with non expressed genes and missing genotype
  genotypes_gene[, is.na(counts_gene[1,])] <- NA
  genotypes_gene[genotypes_gene == -1] <- NA

  ##### Keep genotypes with at least minor_allel_nr number of variants per group; in other case replace them with NAs
  genotypes_gene <- apply(genotypes_gene, 1 ,function(x){
    # x <- genotypes_gene[6,]

    tt <- table(x)

    if( length(tt)==1 )
    return(NULL)
    if( length(tt)==2 ){
      if(any(tt <= minor_allel_nr))
      return(NULL)
      return(x)
      }else{
        if(sum(tt <= minor_allel_nr) >= 2)
        return(NULL)
        x[x == names(tt[tt <= minor_allel_nr])] <- NA
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

names(genotypes_new) <- gene_list
genotypes_new <- MatrixList(genotypes_new)
counts_new <- counts_new[names(genotypes_new)]

data <- new("dmSQTLdata", counts = counts_new, genotypes = genotypes_new, samples = samples)

return(data)

}










