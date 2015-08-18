#' @importClassesFrom IRanges CompressedNumericList
setClassUnion("CompressedNumericListORNULL", c("CompressedNumericList", "NULL"))


##############################################################
#' Object that extends \code{dmSQTLdata} by adding dipsersion.
#' 
#' @slot mean_expression Numeric vector of mean gene expression. When empty equals to NULL.
#' @slot common_dispersion Numeric or NULL when not estimated.
#' @slot tagwise_dispersion CompressedNumericList of estimated tagwise dispersion per gene and genotype or NULL when empty.
setClass("dmSQTLdispersion", 
         contains = "dmSQTLdata",
         representation(mean_expression = "numericORNULL", 
                        common_dispersion = "numericORNULL",
                        tagwise_dispersion = "CompressedNumericListORNULL"), 
         prototype(mean_expression = NULL,
                   common_dispersion = NULL, 
                   tagwise_dispersion = NULL))


##############################################################
setMethod("show", "dmSQTLdispersion", function(object){
  
  callNextMethod(object)
  
  cat("\nSlot \"mean_expression\":\n")
  show_numeric(object@mean_expression)
  
  cat("\nSlot \"common_dispersion\":\n")
  show_numeric(object@common_dispersion)
  
  cat("\nSlot \"tagwise_dispersion\":\n")
  print(object@tagwise_dispersion)
  
})


##############################################################
#' Estimating the dispersion in Dirichlet-multinomial model.
#' 
#' Function that is a wrapper for different optimazation methods used to estimate the dispersion parameters of Dirichlet-multinomial model with maximum likelihood. Parameters that are directly used in the dispersion estimation start with prefix \code{disp_} and the one that are used directly for the proportion estimation start with \code{prop_}.
#' 
#' See details of \code{\link{dmDSdispersion}}
#' 
#' @param x \code{\link{dmSQTLdata}} object with counts and genotypes or any other object that inherits from it.
#' @param ... Parameters needed for the dispersion estimation.
#' @export
setGeneric("dmSQTLdispersion", function(x, ...) standardGeneric("dmSQTLdispersion"))




##############################################################
#' @rdname dmSQTLdispersion
#' @inheritParams dmDSdispersion
#' 
#' @return Returns the \code{\linkS4class{dmSQTLdispersion}} object.
#' @examples 
#' data <- dataSQTL_dmSQTLdata
#' dmSQTLplotData(data)
#' 
#' data <- dmSQTLfilter(data)
#' dmSQTLplotData(data)
#' 
#' \dontrun{
#' ### This part is quite time consuming thus if possible, increase the number of cores.
#' data <- dmSQTLdispersion(data, BPPARAM = BiocParallel::MulticoreParam(workers = 3))
#' }
#' \dontshow{
#' data <- dataSQTL_dmSQTLdispersion
#' }
#' 
#' dmSQTLplotDispersion(data)
#' 
#' @seealso \code{\link[BiocParallel]{bplapply}}
#' @export
setMethod("dmSQTLdispersion", "dmSQTLdata", function(x, mean_expression = TRUE, common_dispersion = FALSE, tagwise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+5), disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "none", disp_prior_df = 10, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  if(mean_expression || (disp_mode == "grid" && disp_moderation == "trended")){
    mean_expression <- dm_estimateMeanExpression(counts = x@counts, BPPARAM = BPPARAM)
  }else{
    mean_expression <- NULL
  }
  
  if(common_dispersion || disp_mode == "grid"){
    common_dispersion <- dmSQTL_estimateCommonDispersion(counts = x@counts, genotypes = x@genotypes, disp_adjust = disp_adjust, disp_interval = disp_interval, disp_tol = 1e+01, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
  }else{
    common_dispersion <- NULL
  }
  
  
  if(tagwise_dispersion){
    
    if(disp_mode == "grid" && !is.null(common_dispersion)){
      cat("!Using common dispersion as disp_init in 'grid' mode!\n")
      disp_init <- common_dispersion
    }
    
    # counts = x@counts; genotypes = x@genotypes
    
    tagwise_dispersion <- dmSQTL_estimateTagwiseDispersion(counts = x@counts, genotypes = x@genotypes, mean_expression = mean_expression, disp_adjust = disp_adjust, disp_mode = disp_mode, disp_interval = disp_interval, disp_tol = disp_tol, disp_init = disp_init, disp_init_weirMoM = disp_init_weirMoM, disp_grid_length = disp_grid_length, disp_grid_range = disp_grid_range, disp_moderation = disp_moderation, disp_prior_df = disp_prior_df, disp_span = disp_span, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
    
  }else{
    tagwise_dispersion <- NULL
  }
  
  
  return(new("dmSQTLdispersion", mean_expression = mean_expression, common_dispersion = common_dispersion, tagwise_dispersion = tagwise_dispersion, counts = x@counts, genotypes = x@genotypes, samples = x@samples))
  
  
})


##############################################################
#' @rdname dmSQTLdispersion
#' @export
setMethod("dmSQTLdispersion", "dmSQTLdispersion", function(x, mean_expression = FALSE, common_dispersion = FALSE, tagwise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = c("optimize", "optim", "constrOptim", "grid")[4], disp_interval = c(0, 1e+5), disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = c("none", "common", "trended")[1], disp_prior_df = 10, disp_span = 0.3, prop_mode = c( "constrOptim", "constrOptimG", "FisherScoring")[2], prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers=1)){
  
  if(mean_expression || (disp_mode == "grid" && disp_moderation == "trended")){
    mean_expression <- dm_estimateMeanExpression(counts = x@counts, BPPARAM = BPPARAM)
  }else{
    mean_expression <- x@mean_expression
  }
  
  if(common_dispersion || disp_mode == "grid"){
    common_dispersion <- dmSQTL_estimateCommonDispersion(counts = x@counts, genotypes = x@genotypes, disp_adjust = disp_adjust, disp_interval = disp_interval, disp_tol = 1e+01, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
  }else{
    common_dispersion <- x@common_dispersion
  }
  
  
  if(tagwise_dispersion){
    
    if(disp_mode == "grid" && !is.null(common_dispersion)){
      cat("!Using common_dispersion =", round(common_dispersion, 2), "as disp_init in 'grid' mode!\n")
      disp_init <- common_dispersion
    }
    
    tagwise_dispersion <- dmSQTL_estimateTagwiseDispersion(counts = x@counts, genotypes = x@genotypes, mean_expression = mean_expression, disp_adjust = disp_adjust, disp_mode = disp_mode, disp_interval = disp_interval, disp_tol = disp_tol, disp_init = disp_init, disp_init_weirMoM = disp_init_weirMoM, disp_grid_length = disp_grid_length, disp_grid_range = disp_grid_range, disp_moderation = disp_moderation, disp_prior_df = disp_prior_df, disp_span = disp_span, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
    
  }else{
    tagwise_dispersion <- x@tagwise_dispersion
  }
  
  
  return(new("dmSQTLdispersion", mean_expression = mean_expression, common_dispersion = common_dispersion, tagwise_dispersion = tagwise_dispersion, counts = x@counts, genotypes = x@genotypes, samples = x@samples))
  
  
})



##############################################################
#' Plot the dispersion - mean trend.
#' 
#' @param x \code{\link{dmSQTLdispersion}} object or any that inherits from it i.e. \code{\link{dmSQTLfit}} or \code{\link{dmSQTLtest}}.
#' @param ... Plotting parameters.
#' @export
setGeneric("dmSQTLplotDispersion", function(x, ...) standardGeneric("dmSQTLplotDispersion"))


##############################################################
#' @rdname dmSQTLplotDispersion
#' @inheritParams dmSQTLplotData
#' @export
setMethod("dmSQTLplotDispersion", "dmSQTLdispersion", function(x, out_dir = NULL){
  
  # tagwise_dispersion = x@tagwise_dispersion; mean_expression = x@mean_expression; nr_features = IRanges::width(x@counts@partitioning); common_dispersion = x@common_dispersion; out_dir = NULL
  
  dmSQTL_plotDispersion(tagwise_dispersion = x@tagwise_dispersion, mean_expression = x@mean_expression, nr_features = IRanges::width(x@counts@partitioning), common_dispersion = x@common_dispersion, out_dir = out_dir)
  
})












































