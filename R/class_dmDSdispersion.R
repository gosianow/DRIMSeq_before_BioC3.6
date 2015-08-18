setClassUnion("numericORNULL", c("numeric", "NULL"))

##############################################################
#' Object that extends \code{dmDSdata} by adding dipsersion.
#' 
#' @slot mean_expression Numeric vector of mean gene expression. When empty equals to NULL.
#' @slot common_dispersion Numeric or NULL when not estimated.
#' @slot tagwise_dispersion Numeric vector of estimated tagwise dispersion.
setClass("dmDSdispersion", 
         contains = "dmDSdata",
         representation(mean_expression = "numericORNULL", 
                        common_dispersion = "numericORNULL",
                        tagwise_dispersion = "numericORNULL"), 
         prototype(mean_expression = NULL,
                   common_dispersion = NULL, 
                   tagwise_dispersion = NULL))


##############################################################
show_numeric <- function(object, nhead = 2, ntail = 2){
  
  if(!is.null(object)){
    
    nl <- length(object)  
    cat(class(object), "of length", length(object), "\n")
    
    if(nl < (nhead + ntail + 1L)) {
      out <- round(object, 2)
    } else {
      dots <- "... ..."
      if(!is.null(names(object)))
      names(dots) <- "... ..."
      out <- c(round(head(object, nhead), 2), dots , round(tail(object, ntail), 2))
    }
    print(out, quote = FALSE, right = TRUE)
    
  }else{
    
    print(object)
    
  }
  
}




setMethod("show", "dmDSdispersion", function(object){
  
  # cat("Slot \"counts\":\n")
  # print(object@counts)
  
  # cat("Slot \"samples\":\n")
  # print(object@samples)
  
  callNextMethod(object)
  
  cat("\nSlot \"mean_expression\":\n")
  show_numeric(object@mean_expression)
  
  cat("\nSlot \"common_dispersion\":\n")
  show_numeric(object@common_dispersion)
  
  cat("\nSlot \"tagwise_dispersion\":\n")
  show_numeric(object@tagwise_dispersion)
  
})


##############################################################
#' Estimating the dispersion in Dirichlet-multinomial model.
#' 
#' Function that is a wrapper for different optimazation methods used to estimate the dispersion parameters of Dirichlet-multinomial model with maximum likelihood. Parameters that are directly used in the dispersion estimation start with prefix \code{disp_} and the one that are used directly for the proportion estimation start with \code{prop_}.
#' 
#' TODO: describe in detail all the methods...
#' 
#' disp_tol: The accuracy defined as \code{tol} in \code{\link{optimize}} when \code{disp_mode = "optimize"}, as \code{factr} in \code{\link{optim}} when \code{disp_mode = "optim"}, as \code{reltol} in \code{\link{constrOptim}} when \code{disp_mode = "constrOptim"}.
#' 
#' @param x \code{\link{dmDSdata}} or \code{\link{dmDSdispersion}} object of counts.
#' @param ... Parameters needed for the dispersion estimation.
#' @export
setGeneric("dmDSdispersion", function(x, ...) standardGeneric("dmDSdispersion"))


##############################################################
#' @rdname dmDSdispersion
#' @param mean_expression Logical. Whether to estimate the mean expression of genes.
#' @param common_dispersion Logical. Whether to estimate the common dispersion.
#' @param tagwise_dispersion Logical. Whether to estimate the tagwise dispersion.
#' @param disp_adjust Logical. Whether to use the adjusted or non-adjusted profile likelihood.
#' @param disp_mode Optimization method used to maximize the profile likelihood. Possible values are \code{"optimize", "optim", "constrOptim", "grid"}. See Details.
#' @param  disp_interval The \code{interval} used by \code{\link{optimize}} function that is used to search for the maximum. Valid when using \code{disp_mode = "optimize"}.
#' @param disp_tol The desired accuracy when estimating dispersion.
#' @param disp_init Initial dispersion.
#' @param disp_init_weirMoM Logical. Whether to use the Weir moment estimator as an initial value for dispersion. If \code{TRUE}, the it overwrites the \code{disp_init} value.
#' @param  disp_grid_length Length of the search grid. Used when \code{disp_mode = "grid"}.
#' @param  disp_grid_range Vector of ends of grid interval.
#' @param disp_moderation Dispersion moderation method. Applies only when \code{disp_mode = "grid"}.
#' @param disp_prior_df Prior degrees of freedom used in moderation.
#' @param disp_span Value from 0 to 1 defining the percentage of genes used in smoothing sliding window when calculating the trend.
#' @param prop_mode Optimization method used to estimate the proportions. Possible values \code{"constrOptim", "constrOptimG"}.
#' @param prop_tol The desired accuracy when estimating dispersion.
#' @param verbose Logical.
#' @param BPPARAM Parallelization method used by \code{\link[BiocParallel]{bplapply}}.
#' 
#' @return Returns the \code{\linkS4class{dmDSdispersion}} object.
#' @examples 
#' data <- dataDS_dmDSdata
#' data <- dmDSfilter(data)
#' \dontrun{
#' data <- dmDSdispersion(data)
#' }
#' \dontshow{
#' data <- dataDS_dmDSdispersion
#' }
#' dmDSplotDispersion(data)
#' 
#' @seealso \code{\link[BiocParallel]{bplapply}}
#' @export
setMethod("dmDSdispersion", "dmDSdata", function(x, mean_expression = TRUE, common_dispersion = TRUE, tagwise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+5), disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "none", disp_prior_df = 10, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers=1)){
  
  if(mean_expression || (disp_mode == "grid" && disp_moderation == "trended")){
    mean_expression <- dm_estimateMeanExpression(counts = x@counts, BPPARAM = BPPARAM)
  }else{
    mean_expression <- NULL
  }
  
  if(common_dispersion || disp_mode == "grid"){
    common_dispersion <- dmDS_estimateCommonDispersion(counts = x@counts, samples = x@samples, disp_adjust = disp_adjust, disp_interval = disp_interval, disp_tol = 1e+01, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
  }else{
    common_dispersion <- NULL
  }
  
  
  if(tagwise_dispersion){
    
    if(disp_mode == "grid" && !is.null(common_dispersion)){
      cat("!Using common dispersion as disp_init in 'grid' mode!\n")
      disp_init <- common_dispersion
    }
    
    tagwise_dispersion <- dmDS_estimateTagwiseDispersion(counts = x@counts, samples = x@samples, mean_expression = mean_expression, disp_adjust = disp_adjust, disp_mode = disp_mode, disp_interval = disp_interval, disp_tol = disp_tol, disp_init = disp_init, disp_init_weirMoM = disp_init_weirMoM, disp_grid_length = disp_grid_length, disp_grid_range = disp_grid_range, disp_moderation = disp_moderation, disp_prior_df = disp_prior_df, disp_span = disp_span, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
    
  }else{
    tagwise_dispersion <- NULL
  }
  
  
  return(new("dmDSdispersion", mean_expression = mean_expression, common_dispersion = common_dispersion, tagwise_dispersion = tagwise_dispersion, counts = x@counts, samples = x@samples))
  
  
})

##############################################################
#' @rdname dmDSdispersion
#' @export
setMethod("dmDSdispersion", "dmDSdispersion", function(x, mean_expression = FALSE, common_dispersion = FALSE, tagwise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+5), disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "none", disp_prior_df = 10, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers=1)){
  
  if(mean_expression || (disp_mode == "grid" && disp_moderation == "trended")){
    mean_expression <- dm_estimateMeanExpression(counts = x@counts, BPPARAM = BPPARAM)
  }else{
    mean_expression <- x@mean_expression
  }
  
  if(common_dispersion || disp_mode == "grid"){
    common_dispersion <- dmDS_estimateCommonDispersion(counts = x@counts, samples = x@samples, disp_adjust = disp_adjust, disp_interval = disp_interval, disp_tol = 1e+01, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
  }else{
    common_dispersion <- x@common_dispersion
  }
  
  
  if(tagwise_dispersion){
    
    if(disp_mode == "grid" && !is.null(common_dispersion)){
      cat("!Using common_dispersion =", round(common_dispersion, 2), "as disp_init in 'grid' mode!\n")
      disp_init <- common_dispersion
    }
    
    tagwise_dispersion <- dmDS_estimateTagwiseDispersion(counts = x@counts, samples = x@samples, mean_expression = mean_expression, disp_adjust = disp_adjust, disp_mode = disp_mode, disp_interval = disp_interval, disp_tol = disp_tol, disp_init = disp_init, disp_init_weirMoM = disp_init_weirMoM, disp_grid_length = disp_grid_length, disp_grid_range = disp_grid_range, disp_moderation = disp_moderation, disp_prior_df = disp_prior_df, disp_span = disp_span, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
    
  }else{
    tagwise_dispersion <- x@tagwise_dispersion
  }
  
  
  return(new("dmDSdispersion", mean_expression = mean_expression, common_dispersion = common_dispersion, tagwise_dispersion = tagwise_dispersion, counts = x@counts, samples = x@samples))
  
  
})

################################################################################
### dmDSplotDispersion
################################################################################

##############################################################
#' Plot the dispersion - mean trend.
#' 
#' @param x \code{\link{dmDSdispersion}} object or any that inherits from it i.e. \code{\link{dmDSfit}} or \code{\link{dmDStest}}.
#' @param ... Plotting parameters.
#' @export
setGeneric("dmDSplotDispersion", function(x, ...) standardGeneric("dmDSplotDispersion"))




##############################################################
#' @rdname dmDSplotDispersion
#' @inheritParams dmDSplotData
#' @examples 
#' dmDSplotDispersion(dataDS_dmDSdispersion)
#' @export
setMethod("dmDSplotDispersion", "dmDSdispersion", function(x, out_dir = NULL){
  
  # tagwise_dispersion = x@tagwise_dispersion; mean_expression = x@mean_expression; nr_features = width(x@counts@partitioning); common_dispersion = x@common_dispersion
  
  dmDS_plotDispersion(tagwise_dispersion = x@tagwise_dispersion, mean_expression = x@mean_expression, nr_features = IRanges::width(x@counts@partitioning), common_dispersion = x@common_dispersion, out_dir = out_dir)
  
})












































