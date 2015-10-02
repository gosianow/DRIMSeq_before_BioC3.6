#' @include class_dmSQTLdata.R
NULL

##############################################################

#' dmSQTLdispersion object
#' 
#' dmSQTLdispersion extends the \code{\linkS4class{dmSQTLdata}} by adding the dispersion estimates of Dirichlet-multinomial distribution used to model the feature (e.g., transcript, exon, exonic bin) ratios for each gene-SNP pair in the sQTL analysis. Result of \code{\link{dmDispersion}}.
#' 
#' @slot mean_expression Numeric vector of mean gene expression.
#' @slot common_dispersion Numeric value of estimated common dispersion.
#' @slot genewise_dispersion List of estimated gene-wise dispersions. Each element of this list is a vector of dispersions estimated for all the genotype blocks assigned to a given gene.
#' @author Malgorzata Nowicka
#' @seealso \code{\link{plotDispersion}}, \code{\linkS4class{dmSQTLdata}}, \code{\linkS4class{dmSQTLfit}}, \code{\linkS4class{dmSQTLtest}}
setClass("dmSQTLdispersion", 
         contains = "dmSQTLdata",
         representation(mean_expression = "numeric", 
                        common_dispersion = "numeric",
                        genewise_dispersion = "list"))



##############################################################

setMethod("show", "dmSQTLdispersion", function(object){
  
  callNextMethod(object)
  
})


##############################################################
# mean_expression = TRUE; common_dispersion = TRUE; genewise_dispersion = TRUE; disp_adjust = TRUE; disp_mode = "grid"; disp_interval = c(0, 1e+4); disp_tol = 1e-08; disp_init = 100; disp_init_weirMoM = TRUE; disp_grid_length = 21; disp_grid_range = c(-10, 10); disp_moderation = "none"; disp_prior_df = 10; disp_span = 0.3; prop_mode = "constrOptimG"; prop_tol = 1e-12; verbose = TRUE; BPPARAM = BiocParallel::MulticoreParam(workers = 10); speed = TRUE


#' @rdname dmDispersion
#' @param speed Logical. If \code{TRUE}, there will be only one dipsersion calculated per gene and it will be assigned to all the blocks matched with this gene. If \code{FALSE}, dispersion is calculated per each gene-block. Such calculation may take a long time, since there can be hundreds of SNPs/blocks per gene.
#' @export
setMethod("dmDispersion", "dmSQTLdata", function(x, mean_expression = TRUE, common_dispersion = TRUE, genewise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+4), disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "none", disp_prior_df = 1, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 1), speed = TRUE){
  
  if(mean_expression || (genewise_dispersion && disp_mode == "grid" && disp_moderation == "trended")){
    mean_expression <- dm_estimateMeanExpression(counts = x@counts, verbose = verbose, BPPARAM = BPPARAM)
  }else{
    mean_expression <- numeric()
  }
  
  
  if(common_dispersion){
    
    ### only one SNP per gene (null model)
    inds <- 1:length(x@genotypes)
    genotypes <- new( "MatrixList", unlistData = matrix(1, nrow = length(x@genotypes), ncol = ncol(x@genotypes)), partitioning = split(inds, factor(names(x@genotypes), levels = names(x@genotypes))) )
    
    common_dispersion <- dmSQTL_estimateCommonDispersion(counts = x@counts, genotypes = genotypes, disp_adjust = disp_adjust, disp_interval = disp_interval, disp_tol = 1e+01, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
    
  }else{
    common_dispersion <- numeric()
  }
  
  
  if(genewise_dispersion){
    
    if(length(common_dispersion)){
      message("! Using common_dispersion = ", round(common_dispersion, 2), " as disp_init !")
      disp_init <- common_dispersion
    }
    
    if(speed){
      
      ### only one SNP per gene (null model)
      inds <- 1:length(x@genotypes)
      genotypes <- new( "MatrixList", unlistData = matrix(1, nrow = length(x@genotypes), ncol = ncol(x@genotypes)), partitioning = split(inds, factor(names(x@genotypes), levels = names(x@genotypes))) )
      
      genewise_dispersion <- dmSQTL_estimateTagwiseDispersion(counts = x@counts, genotypes = genotypes, mean_expression = mean_expression, disp_adjust = disp_adjust, disp_mode = disp_mode, disp_interval = disp_interval, disp_tol = disp_tol, disp_init = disp_init, disp_init_weirMoM = disp_init_weirMoM, disp_grid_length = disp_grid_length, disp_grid_range = disp_grid_range, disp_moderation = disp_moderation, disp_prior_df = disp_prior_df, disp_span = disp_span, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
      
      ### because we keep only one SNP per gene (null model)
      genewise_dispersion <- relist(rep(unlist(genewise_dispersion), times = width(x@genotypes)), x@genotypes@partitioning)
      
    }else{
      
      genewise_dispersion <- dmSQTL_estimateTagwiseDispersion(counts = x@counts, genotypes = x@genotypes, mean_expression = mean_expression, disp_adjust = disp_adjust, disp_mode = disp_mode, disp_interval = disp_interval, disp_tol = disp_tol, disp_init = disp_init, disp_init_weirMoM = disp_init_weirMoM, disp_grid_length = disp_grid_length, disp_grid_range = disp_grid_range, disp_moderation = disp_moderation, disp_prior_df = disp_prior_df, disp_span = disp_span, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
      
    }
    
  }else{
    genewise_dispersion <- list()
  }
  
  
  return(new("dmSQTLdispersion", mean_expression = mean_expression, common_dispersion = common_dispersion, genewise_dispersion = genewise_dispersion, counts = x@counts, genotypes = x@genotypes, blocks = x@blocks, samples = x@samples))
  
  
})


##############################################################

#' @rdname dmDispersion
#' @export
setMethod("dmDispersion", "dmSQTLdispersion", function(x, mean_expression = FALSE, common_dispersion = FALSE, genewise_dispersion = TRUE, disp_adjust = TRUE, disp_mode =  "grid", disp_interval = c(0, 1e+4), disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "none", disp_prior_df = 1, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 1), speed = TRUE){
  
  if(mean_expression || (genewise_dispersion && disp_mode == "grid" && disp_moderation == "trended")){
    mean_expression <- dm_estimateMeanExpression(counts = x@counts, verbose = verbose, BPPARAM = BPPARAM)
  }else{
    mean_expression <- x@mean_expression
  }
  
  
  if(common_dispersion){
    
    ### only one SNP per gene (null model)
    inds <- 1:length(x@genotypes)
    genotypes <- new( "MatrixList", unlistData = matrix(1, nrow = length(x@genotypes), ncol = ncol(x@genotypes)), partitioning = split(inds, factor(names(x@genotypes), levels = names(x@genotypes))) )
    
    common_dispersion <- dmSQTL_estimateCommonDispersion(counts = x@counts, genotypes = genotypes, disp_adjust = disp_adjust, disp_interval = disp_interval, disp_tol = 1e+01, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
    
  }else{
    common_dispersion <- x@common_dispersion
  }
  
  
  if(genewise_dispersion){
    
    if(length(common_dispersion)){
      message("! Using common_dispersion = ", round(common_dispersion, 2), " as disp_init !")
      disp_init <- common_dispersion
    }
    
    if(speed){
      
      ### only one SNP per gene (null model)
      inds <- 1:length(x@genotypes)
      genotypes <- new( "MatrixList", unlistData = matrix(1, nrow = length(x@genotypes), ncol = ncol(x@genotypes)), partitioning = split(inds, factor(names(x@genotypes), levels = names(x@genotypes))) )
      
      genewise_dispersion <- dmSQTL_estimateTagwiseDispersion(counts = x@counts, genotypes = genotypes, mean_expression = mean_expression, disp_adjust = disp_adjust, disp_mode = disp_mode, disp_interval = disp_interval, disp_tol = disp_tol, disp_init = disp_init, disp_init_weirMoM = disp_init_weirMoM, disp_grid_length = disp_grid_length, disp_grid_range = disp_grid_range, disp_moderation = disp_moderation, disp_prior_df = disp_prior_df, disp_span = disp_span, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
      
      ### because we keep only one SNP per gene (null model)
      genewise_dispersion <- relist(rep(unlist(genewise_dispersion), times = width(x@genotypes)), x@genotypes@partitioning)
      
    }else{
      
      genewise_dispersion <- dmSQTL_estimateTagwiseDispersion(counts = x@counts, genotypes = x@genotypes, mean_expression = mean_expression, disp_adjust = disp_adjust, disp_mode = disp_mode, disp_interval = disp_interval, disp_tol = disp_tol, disp_init = disp_init, disp_init_weirMoM = disp_init_weirMoM, disp_grid_length = disp_grid_length, disp_grid_range = disp_grid_range, disp_moderation = disp_moderation, disp_prior_df = disp_prior_df, disp_span = disp_span, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
      
    }
    
  }else{
    genewise_dispersion <- x@genewise_dispersion
  }
  
  
  return(new("dmSQTLdispersion", mean_expression = mean_expression, common_dispersion = common_dispersion, genewise_dispersion = genewise_dispersion, counts = x@counts, genotypes = x@genotypes, blocks = x@blocks, samples = x@samples))
  
  
})


##############################################################

#' @rdname plotDispersion
#' @export
setMethod("plotDispersion", "dmSQTLdispersion", function(x, out_dir = NULL){
  
  dmSQTL_plotDispersion(genewise_dispersion = x@genewise_dispersion, mean_expression = x@mean_expression, nr_features = width(x@counts), common_dispersion = x@common_dispersion, out_dir = out_dir)
  
})










































