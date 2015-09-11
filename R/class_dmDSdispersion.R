#' @include class_dmDSdata.R
NULL

##############################################################

#' Object that extends \code{dmDSdata} by adding dipsersion.
#' 
#' @slot mean_expression Numeric vector of mean gene expression.
#' @slot common_dispersion Numeric value of estimated common dispersion.
#' @slot genewise_dispersion Numeric vector of estimated tagwise dispersion.
setClass("dmDSdispersion", 
         contains = "dmDSdata",
         representation(mean_expression = "numeric", 
                        common_dispersion = "numeric",
                        genewise_dispersion = "numeric"))


##############################################################


#' @rdname dmDSdispersion-class
#' @export
setGeneric("mean_expression", function(x, ...) standardGeneric("mean_expression"))

#' @rdname dmDSdispersion-class
#' @export
setMethod("mean_expression", "dmDSdispersion", function(x){
  
  data.frame(gene_id = names(x@mean_expression), mean_expression = x@mean_expression, stringsAsFactors = FALSE, row.names = NULL)
  
  })


#' @rdname dmDSdispersion-class
#' @export
setGeneric("common_dispersion", function(x, ...) standardGeneric("common_dispersion"))

#' @rdname dmDSdispersion-class
#' @export
setMethod("common_dispersion", "dmDSdispersion", function(x) x@common_dispersion )


#' @rdname dmDSdispersion-class
#' @export
setGeneric("common_dispersion<-", function(x, value) standardGeneric("common_dispersion<-"))

#' @rdname dmDSdispersion-class
#' @export
setMethod("common_dispersion<-", "dmDSdispersion", function(x, value){
  
  return(new("dmDSdispersion", mean_expression = x@mean_expression, common_dispersion = value, genewise_dispersion = x@genewise_dispersion, counts = x@counts, samples = x@samples))
  
  })

#' @rdname dmDSdispersion-class
#' @export
setGeneric("genewise_dispersion", function(x, ...) standardGeneric("genewise_dispersion"))

#' @rdname dmDSdispersion-class
#' @export
setMethod("genewise_dispersion", "dmDSdispersion", function(x){
  
  data.frame(gene_id = names(x@genewise_dispersion), genewise_dispersion = x@genewise_dispersion, stringsAsFactors = FALSE, row.names = NULL)
  
  })


#' @rdname dmDSdispersion-class
#' @export
setGeneric("genewise_dispersion<-", function(x, value) standardGeneric("genewise_dispersion<-"))

#' @rdname dmDSdispersion-class
#' @export
setMethod("genewise_dispersion<-", "dmDSdispersion", function(x, value){
  
  return(new("dmDSdispersion", mean_expression = x@mean_expression, common_dispersion = x@common_dispersion, genewise_dispersion = value, counts = x@counts, samples = x@samples))

  })



##############################################################

setMethod("show", "dmDSdispersion", function(object){
  
  callNextMethod(object)
  
  cat("  mean_expression(), common_dispersion(), genewise_dispersion()\n")
  
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
#' @param x \code{\linkS4class{dmDSdata}} or \code{\linkS4class{dmDSdispersion}} object of counts.
#' @param ... Parameters needed for the dispersion estimation.
#' @export
setGeneric("dmDispersion", function(x, ...) standardGeneric("dmDispersion"))


##############################################################

#' @rdname dmDispersion
#' @param mean_expression Logical. Whether to estimate the mean expression of genes.
#' @param common_dispersion Logical. Whether to estimate the common dispersion.
#' @param genewise_dispersion Logical. Whether to estimate the gene dispersion.
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
#' @return Returns the \code{\linkS4class{dmDSdispersion}} or \code{\linkS4class{dmSQTLdispersion}} object.
#' @seealso \code{\link[BiocParallel]{bplapply}}
#' @examples 
#' d <- dataDS_dmDSdata
#' d <- dmFilter(d)
#' \dontrun{
#' # If possible, increase the number of workers
#' d <- dmDispersion(d, BPPARAM = BiocParallel::MulticoreParam(workers = 1))
#' }
#' @export
setMethod("dmDispersion", "dmDSdata", function(x, mean_expression = TRUE, common_dispersion = TRUE, genewise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+5), disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "none", disp_prior_df = 10, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  if(mean_expression || (genewise_dispersion && disp_mode == "grid" && disp_moderation == "trended")){
    mean_expression <- dm_estimateMeanExpression(counts = x@counts, verbose = verbose, BPPARAM = BPPARAM)
  }else{
    mean_expression <- numeric()
  }
  
  if(common_dispersion){
    common_dispersion <- dmDS_estimateCommonDispersion(counts = x@counts, samples = x@samples, disp_adjust = disp_adjust, disp_interval = disp_interval, disp_tol = 1e+01, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
  }else{
    common_dispersion <- numeric()
  }
  
  
  if(genewise_dispersion){
    
    if(length(common_dispersion)){
      message("! Using common_dispersion = ", round(common_dispersion, 2), " as disp_init !")
      disp_init <- common_dispersion
    }
    
    genewise_dispersion <- dmDS_estimateTagwiseDispersion(counts = x@counts, samples = x@samples, mean_expression = mean_expression, disp_adjust = disp_adjust, disp_mode = disp_mode, disp_interval = disp_interval, disp_tol = disp_tol, disp_init = disp_init, disp_init_weirMoM = disp_init_weirMoM, disp_grid_length = disp_grid_length, disp_grid_range = disp_grid_range, disp_moderation = disp_moderation, disp_prior_df = disp_prior_df, disp_span = disp_span, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
    
  }else{
    genewise_dispersion <- numeric()
  }
  
  
  return(new("dmDSdispersion", mean_expression = mean_expression, common_dispersion = common_dispersion, genewise_dispersion = genewise_dispersion, counts = x@counts, samples = x@samples))
  
  
})

##############################################################
#' @rdname dmDispersion
#' @export
setMethod("dmDispersion", "dmDSdispersion", function(x, mean_expression = FALSE, common_dispersion = FALSE, genewise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+5), disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "none", disp_prior_df = 10, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  if(mean_expression || (genewise_dispersion && disp_mode == "grid" && disp_moderation == "trended")){
    mean_expression <- dm_estimateMeanExpression(counts = x@counts, verbose = verbose, BPPARAM = BPPARAM)
  }else{
    mean_expression <- x@mean_expression
  }
  
  if(common_dispersion){
    common_dispersion <- dmDS_estimateCommonDispersion(counts = x@counts, samples = x@samples, disp_adjust = disp_adjust, disp_interval = disp_interval, disp_tol = 1e+01, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
  }else{
    common_dispersion <- x@common_dispersion
  }
  
  
  if(genewise_dispersion){
    
    if(length(common_dispersion)){
      message("! Using common_dispersion = ", round(common_dispersion, 2), " as disp_init !")
      disp_init <- common_dispersion
    }
    
    genewise_dispersion <- dmDS_estimateTagwiseDispersion(counts = x@counts, samples = x@samples, mean_expression = mean_expression, disp_adjust = disp_adjust, disp_mode = disp_mode, disp_interval = disp_interval, disp_tol = disp_tol, disp_init = disp_init, disp_init_weirMoM = disp_init_weirMoM, disp_grid_length = disp_grid_length, disp_grid_range = disp_grid_range, disp_moderation = disp_moderation, disp_prior_df = disp_prior_df, disp_span = disp_span, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
    
  }else{
    genewise_dispersion <- x@genewise_dispersion
  }
  
  
  return(new("dmDSdispersion", mean_expression = mean_expression, common_dispersion = common_dispersion, genewise_dispersion = genewise_dispersion, counts = x@counts, samples = x@samples))
  
  
})


##############################################################

#' Plot the dispersion - mean trend.
#' 
#' @param x \code{\linkS4class{dmDSdispersion}} or \code{\linkS4class{dmSQTLdispersion}} object or any that inherits from it.
#' @param ... Plotting parameters.
#' @export
setGeneric("plotDispersion", function(x, ...) standardGeneric("plotDispersion"))




##############################################################

#' @rdname plotDispersion
#' @inheritParams plotData
#' @examples 
#' d <- dataDS_dmDSdispersion
#' plotDispersion(d)
#' plot(d)
#' @export
setMethod("plotDispersion", "dmDSdispersion", function(x, out_dir = NULL){
  
  if(!length(x@genewise_dispersion) == length(x@counts))
  stop("Genewise dispersion must be estimated for each gene!")
  if(!length(x@genewise_dispersion) == length(x@mean_expression))
  stop("Mean expression must be estimated for each gene!")
  
  dmDS_plotDispersion(genewise_dispersion = x@genewise_dispersion, mean_expression = x@mean_expression, nr_features = width(x@counts), common_dispersion = x@common_dispersion, out_dir = out_dir)
  
})


##############################################################

#' @rdname dmDSdispersion-class
#' @export
setMethod("plot", "dmDSdispersion", function(x, out_dir = NULL){
  
  plotDispersion(x, out_dir = out_dir)
  
})









































