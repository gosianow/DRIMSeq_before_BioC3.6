################################################################################
#' Object that extends \code{dmDSfit} by adding the test results.
#' 
#' @slot table data.frame with \code{gene_id} - gene IDs, \code{lr} - likelihood ratio statistics, \code{df} - degrees of freedom, \code{pvalue} - p-values and \code{adj_pvalue} - Benjamini & Hochberg adjusted p-values.
setClass("dmDStest", 
         contains = "dmDSfit",
         representation(table = "data.frame"))

################################################################################
setMethod("show", "dmDStest", function(object){
  
  callNextMethod(object)
  
  cat("\nSlot \"table\":\n")
  show_matrix(object@table)
  
  
})

################################################################################
#' Likelihood ratio test between full and null model.
#' 
#' @param x \code{\link{dmDStest}} object.
#' @param ... Parameters needed for the likelihood ratio test.
#' @export
setGeneric("dmDStest", function(x, ...) standardGeneric("dmDStest"))

################################################################################
#' @rdname dmDStest
#' @return This function returns a \code{\link{dmDStest}} object with an additional slot \code{table} which is sorted by significance and contains  \code{gene_id} - gene IDs, \code{lr} - likelihood ratio statistics, \code{df} - degrees of freedom, \code{pvalue} - p-values and \code{adj_pvalue} - Benjamini & Hochberg adjusted p-values.
#' @examples 
#' data <- dataDS_dmDSdispersion
#' data <- dmDSfit(data)
#' data <- dmDStest(data)
#' dmDSplotTest(data)
#' @export
setMethod("dmDStest", "dmDSfit", function(x){
  
  
  table <- dmDS_test(stats_full = x@fit_full@metadata, stats_null = x@fit_null@metadata)
  
  
  return(new("dmDStest", table = table, dispersion = x@dispersion, fit_full = x@fit_full, fit_null = x@fit_null, mean_expression = x@mean_expression, common_dispersion = x@common_dispersion, tagwise_dispersion = x@tagwise_dispersion, counts = x@counts, samples = x@samples))
  
  
})


################################################################################
#' Plot the histogram of p-values.
#' 
#' @param x \code{\link{dmDStest}} object.
#' @param ... Plotting parameters.
#' @export
setGeneric("dmDSplotTest", function(x, ...) standardGeneric("dmDSplotTest"))



################################################################################
#' @rdname dmDSplotTest
#' @inheritParams dmDSplotFit
#' @export
setMethod("dmDSplotTest", "dmDStest", function(x, out_dir = NULL){
  
  dm_plotTable(x@table, out_dir = out_dir)
  
})

################################################################################
#' @rdname dmDSplotFit
#' @export
setMethod("dmDSplotFit", "dmDStest", function(x, gene_id, plot_type = c("barplot", "boxplot1", "boxplot2", "lineplot", "ribbonplot")[1], order = TRUE, plot_full = TRUE, plot_null = TRUE, out_dir = NULL){
  
  stopifnot(plot_type %in% c("barplot", "boxplot1", "boxplot2", "lineplot", "ribbonplot"))
  
  dmDS_plotFit(gene_id = gene_id, counts = x@counts, samples = x@samples, dispersion = slot(x, x@dispersion), proportions_full = x@fit_full, proportions_null = x@fit_full, table = x@table, plot_type = plot_type, order = order, plot_full = plot_full, plot_null = plot_null, out_dir = out_dir)
  
  
})













































