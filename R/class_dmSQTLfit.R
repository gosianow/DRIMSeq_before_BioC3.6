#' @include class_dmSQTLdispersion.R
NULL

################################################################################

#' Object that extends \code{dmSQTLdispersion} by adding fitting.
#' 
#' @slot dispersion Character specifying which type of dispersion was used for fitting.
#' @slot fit_full List of \code{\linkS4class{MatrixList}} objects. Each of them contains the per gene fitting results that include the group proportions (groups are defined by genotypes) of feature expression, full model likelihoods and degrees of freedom.
setClass("dmSQTLfit", 
  contains = "dmSQTLdispersion",
  representation(dispersion = "character",
    fit_full = "list"))

################################################################################

setMethod("show", "dmSQTLfit", function(object){
  
  callNextMethod(object)
  
  cat("\nSlot \"dispersion\":\n")
  print(object@dispersion)
  
  cat("\nSlot \"fit_full\":\n")
  show_MatrixList_list(object@fit_full)
  
  
  })


################################################################################

#' @rdname dmFit
#' @export
setMethod("dmFit", "dmSQTLdispersion", function(x, dispersion = "genewise_dispersion", prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
   
  fit_full <- dmSQTL_fitOneModel(counts = x@counts, genotypes = x@genotypes, dispersion = slot(x, dispersion), model = "full", prop_mode = prop_mode, prop_tol = prop_tol, verbose = TRUE, BPPARAM = BPPARAM)

 
  return(new("dmSQTLfit", dispersion = dispersion, fit_full = fit_full, mean_expression = x@mean_expression, common_dispersion = x@common_dispersion, genewise_dispersion = x@genewise_dispersion, counts = x@counts, genotypes = x@genotypes, samples = x@samples))
  
  
  })



################################################################################

#' @rdname plotFit
#' @export
setMethod("plotFit", "dmSQTLfit", function(x, gene_id, snp_id, plot_type = "boxplot1", order = TRUE, plot_full = TRUE, out_dir = NULL){
  
  stopifnot(plot_type %in% c("barplot", "boxplot1", "boxplot2", "lineplot", "ribbonplot"))
  
  dmSQTL_plotFit(gene_id = gene_id, snp_id = snp_id, counts = x@counts, genotypes = x@genotypes, samples = x@samples, dispersion = slot(x, x@dispersion), fit_full = x@fit_full, fit_null = NULL, table = NULL, plot_type = plot_type, order = order, plot_full = plot_full, plot_null = FALSE, out_dir = out_dir)
  
  
  })




##############################################################

#' @rdname dmSQTLfit-class
#' @export
setMethod("plot", "dmSQTLfit", function(x, gene_id, snp_id, plot_type = "boxplot1", order = TRUE, plot_full = TRUE, out_dir = NULL){
  
  plotFit(x, gene_id, snp_id, plot_type = plot_type, order = order, plot_full = plot_full, out_dir = out_dir)
  
})













































