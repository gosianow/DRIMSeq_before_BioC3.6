#' @include class_dmDSfit.R
NULL

################################################################################

#' Object that extends \code{dmDSfit} by adding the test results.
#' 
#' @slot table data.frame with \code{gene_id} - gene IDs, \code{lr} - likelihood ratio statistics, \code{df} - degrees of freedom, \code{pvalue} - p-values and \code{adj_pvalue} - Benjamini & Hochberg adjusted p-values.
setClass("dmDSLRT", 
         contains = "dmDSfit",
         representation(table = "data.frame"))

################################################################################

setMethod("show", "dmDSLRT", function(object){
  
  callNextMethod(object)
  
  cat("\nSlot \"table\":\n")
  show_matrix(object@table, nhead = 4, ntail = 4)
  
  
})

################################################################################

#' Likelihood ratio test between full and null model.
#' 
#' @param x \code{\linkS4class{dmDSfit}} or \code{\linkS4class{dmSQTLfit}} object.
#' @param ... Parameters needed for the likelihood ratio test.
#' @export
setGeneric("dmLRT", function(x, ...) standardGeneric("dmLRT"))

################################################################################

#' @rdname dmLRT
#' @return This function returns a \code{\linkS4class{dmDSLRT}} or \code{\linkS4class{dmSQTLLRT}} object with an additional slot \code{table} which is sorted by significance and contains  \code{gene_id} - gene IDs, \code{lr} - likelihood ratio statistics, \code{df} - degrees of freedom, \code{pvalue} - p-values and \code{adj_pvalue} - Benjamini & Hochberg adjusted p-values.
#' @export
setMethod("dmLRT", "dmDSfit", function(x){
  
  
  table <- dmDS_test(stats_full = x@fit_full@metadata, stats_null = x@fit_null@metadata)
  
  
  return(new("dmDSLRT", table = table, dispersion = x@dispersion, fit_full = x@fit_full, fit_null = x@fit_null, mean_expression = x@mean_expression, common_dispersion = x@common_dispersion, genewise_dispersion = x@genewise_dispersion, counts = x@counts, samples = x@samples))
  
  
})


################################################################################

#' Plot the histogram of p-values.
#' 
#' @param x \code{\linkS4class{dmDSLRT}} or \code{\linkS4class{dmSQTLLRT}} object.
#' @param ... Plotting parameters.
#' @export
setGeneric("plotLRT", function(x, ...) standardGeneric("plotLRT"))



################################################################################

#' @rdname plotLRT
#' @inheritParams plotFit
#' @export
setMethod("plotLRT", "dmDSLRT", function(x, out_dir = NULL){
  
  dm_plotTable(x@table, out_dir = out_dir)
  
})

################################################################################

#' @rdname plotFit
#' @export
setMethod("plotFit", "dmDSLRT", function(x, gene_id, plot_type = c("barplot", "boxplot1", "boxplot2", "lineplot", "ribbonplot")[1], order = TRUE, plot_full = TRUE, plot_null = TRUE, out_dir = NULL){
  
  stopifnot(plot_type %in% c("barplot", "boxplot1", "boxplot2", "lineplot", "ribbonplot"))
  
  dmDS_plotFit(gene_id = gene_id, counts = x@counts, samples = x@samples, dispersion = slot(x, x@dispersion), proportions_full = x@fit_full, proportions_null = x@fit_full, table = x@table, plot_type = plot_type, order = order, plot_full = plot_full, plot_null = plot_null, out_dir = out_dir)
  
  
})













































