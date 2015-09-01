#' @include class_dmDSdispersion.R
NULL

################################################################################

#' Object that extends \code{dmDSdispersion} class by adding fitting.
#' 
#' @slot dispersion Character specifying which type of dispersion was used for fitting: common_dispersion or genewise_dispersion.
#' @slot fit_full \code{\linkS4class{MatrixList}} object containing the per group proportions of feature expression and likelihoods.
setClass("dmDSfit", 
         contains = "dmDSdispersion",
         representation(dispersion = "character",
                        fit_full = "MatrixList"))


################################################################################


#' @rdname dmDSfit-class
#' @export
setGeneric("proportions", function(x, ...) standardGeneric("proportions"))

#' @rdname dmDSfit-class
#' @export
setMethod("proportions", "dmDSfit", function(x){
  
  data.frame(gene_id = rep(names(x@counts), width(x@counts)), feature_id = rownames(x@counts@unlistData), x@fit_full@unlistData, stringsAsFactors = FALSE, row.names = NULL)
  
    })


#' @rdname dmDSfit-class
#' @export
setGeneric("statistics", function(x, ...) standardGeneric("statistics"))

#' @rdname dmDSfit-class
#' @export
setMethod("statistics", "dmDSfit", function(x){
  
  df <- data.frame(gene_id = names(x@counts), x@fit_full@metadata, stringsAsFactors = FALSE, row.names = NULL)
  colnames(df)[-1] <- paste0("lik_", colnames(x@fit_full@metadata))
  return(df)
  
    })



################################################################################

setMethod("show", "dmDSfit", function(object){
  
  callNextMethod(object)
  
  cat("  proportions(), statistics()\n")
  
  
})

################################################################################

#' Estimating the proportions and likelihoods of Dirichlet-multinomial.
#' 
#' @param x \code{\linkS4class{dmDSdispersion}} or \code{\linkS4class{dmSQTLdispersion}} object or any that inherits from them.
#' @param ... Parameters needed for the proportion estimation.
#' @export
setGeneric("dmFit", function(x, ...) standardGeneric("dmFit"))


################################################################################

#' @rdname dmFit
#' @inheritParams dmDispersion
#' @param dispersion Characted defining which dispersion should be used for fitting. Possible values \code{"genewise_dispersion"}, \code{"common_dispersion"}
#' @examples 
#' d <- dataDS_dmDSdispersion
#' d <- dmFit(d)
#' @export
setMethod("dmFit", "dmDSdispersion", function(x, dispersion = "genewise_dispersion", prop_mode = c("constrOptim", "constrOptimG", "FisherScoring")[2], prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  
  fit_full <- dmDS_fitOneModel(counts = x@counts, samples = x@samples, dispersion = slot(x, dispersion), model = "full", prop_mode = prop_mode, prop_tol = prop_tol, verbose = TRUE, BPPARAM = BPPARAM)
  
  
  return(new("dmDSfit", dispersion = dispersion, fit_full = fit_full, mean_expression = x@mean_expression, common_dispersion = x@common_dispersion, genewise_dispersion = x@genewise_dispersion, counts = x@counts, samples = x@samples))
  
  
})


################################################################################

#' Plot the estimated proportions.
#' 
#' @param x \code{\linkS4class{dmDSfit}} or \code{\linkS4class{dmSQTLfit}} object or any that inherits from them.
#' @param ... Plotting parameters.
#' @export
setGeneric("plotFit", function(x, ...) standardGeneric("plotFit"))


################################################################################

#' @rdname plotFit
#' @inheritParams plotDispersion
#' @param gene_id Vector of gene IDs to be plotted.
#' @param plot_type Character defining type of the plot produced. Possible values \code{"barplot"}, \code{"boxplot1"}, \code{"boxplot2"}, \code{"lineplot"}, \code{"ribbonplot"}.
#' @param order Logical. Whether to plot the features ordered by their expression.
#' @param plot_full Logical. Whether to plot the proportions estimated by the full model.
#' @param plot_main Logical. Whether to plot the plot title.
#' @examples 
#' d <- dataDS_dmDSdispersion
#' 
#' # If possible, increase the number of workers
#' d <- dmFit(d, BPPARAM = BiocParallel::MulticoreParam(workers = 1))
#' 
#' gene_id <- names(d)[1]
#' 
#' plotFit(d, gene_id = gene_id)
#' plot(d, gene_id = gene_id, plot_type = "lineplot", plot_full = FALSE)
#' 
#' d <- dmLRT(d)
#' 
#' results <- results(d)
#' 
#' gene_id <- results$gene_id[1:3]
#' 
#' plotFit(d, gene_id = gene_id)
#' 
#' 
#' @export
setMethod("plotFit", "dmDSfit", function(x, gene_id, plot_type = "barplot", order = TRUE, plot_full = TRUE, plot_main = TRUE, out_dir = NULL){
  
  stopifnot(plot_type %in% c("barplot", "boxplot1", "boxplot2", "lineplot", "ribbonplot"))
  
  dmDS_plotFit(gene_id = gene_id, counts = x@counts, samples = x@samples, dispersion = slot(x, x@dispersion), proportions_full = x@fit_full, proportions_null = NULL, table = NULL, plot_type = plot_type, order = order, plot_full = plot_full, plot_null = FALSE, plot_main = plot_main, out_dir = out_dir)
  
  
})



##############################################################

#' @rdname dmDSfit-class
#' @examples 
#' d <- dataDS_dmDSdispersion
#' 
#' # If possible, increase the number of workers
#' d <- dmFit(d, BPPARAM = BiocParallel::MulticoreParam(workers = 1))
#' 
#' gene_id <- names(d)[1]
#' 
#' plotFit(d, gene_id = gene_id)
#' plot(d, gene_id = gene_id, plot_type = "lineplot", plot_full = FALSE)
#' 
#' @export
setMethod("plot", "dmDSfit", function(x, gene_id, plot_type = "barplot", order = TRUE, plot_full = TRUE, plot_main = TRUE, out_dir = NULL){
  
  plotFit(x, gene_id, plot_type = plot_type, order = order, plot_full = plot_full, plot_main = plot_main, out_dir = out_dir)
  
})















































