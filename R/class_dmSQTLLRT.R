#' @include class_dmSQTLfit.R
NULL

################################################################################

#' Object that extends \code{dmSQTLfit} by adding the test results.
#' 
#' @slot table data.frame with \code{gene_id} - gene IDs, \code{snp_id} - SNP IDs, \code{lr} - likelihood ratio statistics, \code{df} - degrees of freedom, \code{pvalue} - p-values and \code{adj_pvalue} - Benjamini & Hochberg adjusted p-values.
setClass("dmSQTLLRT", 
         contains = "dmSQTLfit",
         representation(table = "data.frame"))

################################################################################

setMethod("show", "dmSQTLLRT", function(object){
  
  callNextMethod(object)
  
  cat("\nSlot \"table\":\n")
  show_matrix(object@table, nhead = 4, ntail = 4)
  
  
})


################################################################################

#' @rdname dmLRT
#' @export
setMethod("dmLRT", "dmSQTLfit", function(x, BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  table <- dmSQTL_test(fit_full = x@fit_full, fit_null = x@fit_null, BPPARAM = BPPARAM)
  
  return(new("dmSQTLLRT", table = table, dispersion = x@dispersion, fit_full = x@fit_full, fit_null = x@fit_null, mean_expression = x@mean_expression, common_dispersion = x@common_dispersion, tagwise_dispersion = x@tagwise_dispersion, counts = x@counts, genotypes = x@genotypes, samples = x@samples))
  
  
})


################################################################################

#' @rdname plotLRT
#' @export
setMethod("plotLRT", "dmSQTLLRT", function(x, out_dir = NULL){
  
  dm_plotTable(x@table, out_dir = out_dir)
  
})


################################################################################

#' @rdname plotFit
#' @export
setMethod("plotFit", "dmSQTLLRT", function(x, gene_id, snp_id, plot_type = "boxplot1", order = TRUE, plot_full = TRUE, plot_null = TRUE, out_dir = NULL){
  
  stopifnot(plot_type %in% c("barplot", "boxplot1", "boxplot2", "lineplot", "ribbonplot"))
  
  dmSQTL_plotFit(gene_id = gene_id, snp_id = snp_id, counts = x@counts, genotypes = x@genotypes, samples = x@samples, dispersion = slot(x, x@dispersion), fit_full = x@fit_full, fit_null = x@fit_null, table = x@table, plot_type = plot_type, order = order, plot_full = plot_full, plot_null = plot_null, out_dir = out_dir)
  
  
})













































