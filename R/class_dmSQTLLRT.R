#' @include class_dmSQTLfit.R
NULL

################################################################################

#' Object that extends \code{dmSQTLfit} by adding the test results.
#' 
#' @slot fit_null list
#' @slot results data.frame with \code{gene_id} - gene IDs, \code{snp_id} - SNP IDs, \code{lr} - likelihood ratio statistics, \code{df} - degrees of freedom, \code{pvalue} - p-values and \code{adj_pvalue} - Benjamini & Hochberg adjusted p-values.
setClass("dmSQTLLRT", 
         contains = "dmSQTLfit",
         representation(fit_null = "list",
          results = "data.frame"))

################################################################################

setMethod("show", "dmSQTLLRT", function(object){
  
  callNextMethod(object)
  
  cat("\nSlot \"fit_null\":\n")
  show_MatrixList_list(object@fit_null)
  
  cat("\nSlot \"results\":\n")
  show_matrix(object@results, nhead = 5, ntail = 5)
  
  
})


################################################################################

#' @rdname dmLRT
#' @export
setMethod("dmLRT", "dmSQTLfit", function(x, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  fit_null <- dmSQTL_fitOneModel(counts = x@counts, genotypes = x@genotypes, dispersion = slot(x, x@dispersion), model = "null", prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
  
  results <- dmSQTL_test(fit_full = x@fit_full, fit_null = fit_null, BPPARAM = BPPARAM)
  
  return(new("dmSQTLLRT", fit_null = fit_null, results = results, dispersion = x@dispersion, fit_full = x@fit_full, mean_expression = x@mean_expression, common_dispersion = x@common_dispersion, genewise_dispersion = x@genewise_dispersion, counts = x@counts, genotypes = x@genotypes, samples = x@samples))
  
  
})


################################################################################

#' @rdname plotLRT
#' @export
setMethod("plotLRT", "dmSQTLLRT", function(x, out_dir = NULL){
  
  dm_plotTable(x@results, out_dir = out_dir)
  
})


################################################################################

#' @rdname plotFit
#' @export
setMethod("plotFit", "dmSQTLLRT", function(x, gene_id, snp_id, plot_type = "boxplot1", order = TRUE, plot_full = TRUE, plot_null = TRUE, out_dir = NULL){
  
  stopifnot(plot_type %in% c("barplot", "boxplot1", "boxplot2", "lineplot", "ribbonplot"))
  
  dmSQTL_plotFit(gene_id = gene_id, snp_id = snp_id, counts = x@counts, genotypes = x@genotypes, samples = x@samples, dispersion = slot(x, x@dispersion), fit_full = x@fit_full, fit_null = x@fit_null, table = x@results, plot_type = plot_type, order = order, plot_full = plot_full, plot_null = plot_null, out_dir = out_dir)
  
  
})













































