#' @include class_dmFit.R
NULL

################################################################################
### dmDSfit
################################################################################


#' Object that extends \code{dmDSdispersion} by adding fitting.
#' 
#' @slot dispersion Character specifying which type of dispersion was used for fitting.
#' @slot fit_full \code{\linkS4class{dmFit}} object containing the per group proportions of feature expression and full model likelihoods and degrees of freedom.
#' @slot fit_null \code{\linkS4class{dmFit}} object containing the pooled proportions of feature expression and null model likelihoods and degrees of freedom.
setClass("dmDSfit", 
         contains = "dmDSdispersion",
         representation(dispersion = "character",
                        fit_full = "dmFit",
                        fit_null = "dmFit"))


setMethod("show", "dmDSfit", function(object){
  
  callNextMethod(object)
  
  cat("\nSlot \"dispersion\":\n")
  print(object@dispersion)
  
  cat("\nSlot \"fit_full\":\n")
  print(object@fit_full)
  
  cat("\nSlot \"fit_null\":\n")
  print(object@fit_null)
  
  
})


#' Estimating the proportions and likelihoods of Dirichlet-multinomial full and null models.
#' 
#' @param x \code{\link{dmDSdispersion}} object or any that inherits from it i.e. \code{\link{dmDSfit}} or \code{\link{dmDStest}}.
#' @param ... Parameters needed for the proportion estimation.
#' @export
setGeneric("dmDSfit", function(x, ...) standardGeneric("dmDSfit"))



#' @rdname dmDSfit
#' @inheritParams dmDSdispersion
#' @param dispersion Characted defining which dispersion should be used for fitting. Possible values \code{"tagwise_dispersion", "common_dispersion"}
#' @examples 
#' data <- dataDS_dmDSdispersion
#' data <- dmDSfit(data)
#' dmDSplotFit(data, gene_id = "FBgn0001316", plot_type = "barplot")
#' @export
setMethod("dmDSfit", "dmDSdispersion", function(x, dispersion = "tagwise_dispersion", prop_mode = c("constrOptim", "constrOptimG", "FisherScoring")[2], prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers=1)){
  
  
  fit_full <- dmDS_fitOneModel(counts = x@counts, samples = x@samples, dispersion = slot(x, dispersion), model = "full", prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
  
  
  fit_null <- dmDS_fitOneModel(counts = x@counts, samples = x@samples, dispersion = slot(x, dispersion), model = "null", prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
  
  
  return(new("dmDSfit", dispersion = dispersion, fit_full = fit_full, fit_null = fit_null, mean_expression = x@mean_expression, common_dispersion = x@common_dispersion, tagwise_dispersion = x@tagwise_dispersion, counts = x@counts, samples = x@samples))
  
  
})


################################################################################
### dmDSplotFit
################################################################################

#' Plot the estimated proportions.
#' 
#' @param x \code{\link{dmDSfit}} object or any that inherits from it i.e. \code{\link{dmDStest}}.
#' @param ... Plotting parameters.
#' @export
setGeneric("dmDSplotFit", function(x, ...) standardGeneric("dmDSplotFit"))



#' @rdname dmDSplotFit
#' @inheritParams dmDSplotDispersion
#' @param gene_id Vector of gene IDs to be plotted.
#' @param plot_type Character defining type of the plot produced. Possible values \code{"barplot", "boxplot1", "boxplot2", "lineplot", "ribbonplot"}.
#' @param order Logical. Whether to plot the features ordered by their expression.
#' @param plot_full Logical. Whether to plot the proportions estimated by the full model.
#' @param plot_null Logical. Whether to plot the proportions estimated by the null model.
#' @examples 
#' data <- dataDS_dmDSdispersion
#' data <- dmDSfit(data)
#' dmDSplotFit(data, gene_id = "FBgn0001316", plot_type = "barplot")
#' dmDSplotFit(data, gene_id = "FBgn0001316", plot_type = "barplot", plot_full = FALSE, plot_null = FALSE)
#' @export
setMethod("dmDSplotFit", "dmDSfit", function(x, gene_id, plot_type = "barplot", order = TRUE, plot_full = TRUE, plot_null = TRUE, out_dir = NULL){
  
  
  dmDS_plotFit(gene_id = gene_id, counts = x@counts, samples = x@samples, dispersion = slot(x, x@dispersion), proportions_full = x@fit_full@proportions, proportions_null = x@fit_null@proportions, table = NULL, plot_type = plot_type, order = order, plot_full = plot_full, plot_null = plot_null, out_dir = out_dir)
  
  
})

















































