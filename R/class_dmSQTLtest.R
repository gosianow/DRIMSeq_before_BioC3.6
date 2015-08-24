################################################################################

#' Object that extends \code{dmSQTLfit} by adding the test results.
#' 
#' @slot table data.frame with \code{gene_id} - gene IDs, \code{snp_id} - SNP IDs, \code{lr} - likelihood ratio statistics, \code{df} - degrees of freedom, \code{pvalue} - p-values and \code{adj_pvalue} - Benjamini & Hochberg adjusted p-values.
setClass("dmSQTLtest", 
         contains = "dmSQTLfit",
         representation(table = "data.frame"))

################################################################################
setMethod("show", "dmSQTLtest", function(object){
  
  callNextMethod(object)
  
  cat("\nSlot \"table\":\n")
  show_matrix(object@table)
  
  
})


################################################################################
#' Likelihood ratio test between full and null model.
#' 
#' @param x \code{\link{dmSQTLtest}} object.
#' @param ... Parameters needed for the likelihood ratio test.
#' @export
setGeneric("dmSQTLtest", function(x, ...) standardGeneric("dmSQTLtest"))



################################################################################
#' @rdname dmSQTLtest
#' @inheritParams dmSQTLfit
#' @return This function returns a \code{\link{dmSQTLtest}} object with an additional slot \code{table} which is sorted by significance and contains  \code{gene_id} - gene IDs, \code{snp_id} - SNP IDs, \code{lr} - likelihood ratio statistics, \code{df} - degrees of freedom, \code{pvalue} - p-values and \code{adj_pvalue} - Benjamini & Hochberg adjusted p-values.
#' @examples 
#' data <- dataSQTL_dmSQTLdispersion
#' data <- dmSQTLfit(data)
#'
#' data <- dmSQTLtest(data)
#' dmSQTLplotTest(data)
#' 
#' snp_id <- "snp_19_48981946"
#' gene_id <- "ENSG00000105443.8"
#' dmSQTLplotFit(data, gene_id, snp_id)
#' 
#' @export
setMethod("dmSQTLtest", "dmSQTLfit", function(x, BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  # fit_full = x@fit_full; fit_null = x@fit_null
  
  table <- dmSQTL_test(fit_full = x@fit_full, fit_null = x@fit_null, BPPARAM = BPPARAM)
  
  
  return(new("dmSQTLtest", table = table, dispersion = x@dispersion, fit_full = x@fit_full, fit_null = x@fit_null, mean_expression = x@mean_expression, common_dispersion = x@common_dispersion, tagwise_dispersion = x@tagwise_dispersion, counts = x@counts, genotypes = x@genotypes, samples = x@samples))
  
  
})


################################################################################
#' Plot the histogram of p-values.
#' 
#' @param x \code{\link{dmSQTLtest}} object.
#' @param ... Plotting parameters.
#' @export
setGeneric("dmSQTLplotTest", function(x, ...) standardGeneric("dmSQTLplotTest"))



################################################################################
#' @rdname dmSQTLplotTest
#' @inheritParams dmSQTLplotFit
#' @export
setMethod("dmSQTLplotTest", "dmSQTLtest", function(x, out_dir = NULL){
  
  dm_plotTable(x@table, out_dir = out_dir)
  
})


################################################################################
#' @rdname dmSQTLplotFit
#' @export
setMethod("dmSQTLplotFit", "dmSQTLtest", function(x, gene_id, snp_id, plot_type = "boxplot1", order = TRUE, plot_full = TRUE, plot_null = TRUE, out_dir = NULL){
  
  stopifnot(plot_type %in% c("barplot", "boxplot1", "boxplot2", "lineplot", "ribbonplot"))
  
  dmSQTL_plotFit(gene_id = gene_id, snp_id = snp_id, counts = x@counts, genotypes = x@genotypes, samples = x@samples, dispersion = slot(x, x@dispersion), fit_full = x@fit_full, fit_null = x@fit_null, table = x@table, plot_type = plot_type, order = order, plot_full = plot_full, plot_null = plot_null, out_dir = out_dir)
  
  
})













































