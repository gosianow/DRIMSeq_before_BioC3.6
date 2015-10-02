#' @include class_dmSQTLdispersion.R class_dmDSfit.R
NULL

################################################################################

#' dmSQTLfit object
#' 
#' dmSQTLfit extends the \code{\linkS4class{dmDSdispersion}} class by adding the full model Dirichlet-multinomial feature proportion estimates needed for the sQTL analysis. Feature ratios are estimated for each gene and each group that is defined by different SNPs/blocks. Result of dmFit.
#' 
#' @slot dispersion Character specifying which type of dispersion was used for fitting: \code{"common_dispersion"} or \code{"genewise_dispersion"}.
#' @slot fit_full List of \code{\linkS4class{MatrixList}} objects. Each element of this list contains the full model proportion estimates for all the blocks associated with a given gene. Columns of MatrixLists correspond to 3 genotypes (0,1,2). The full model likelihoods are stored in \code{metadata} slot.
#' @author Malgorzata Nowicka
#' @seealso \code{\link{plotFit}}, \code{\linkS4class{dmSQTLdata}}, \code{\linkS4class{dmSQTLdispersion}}, \code{\linkS4class{dmSQTLtest}}
setClass("dmSQTLfit", 
         contains = "dmSQTLdispersion",
         representation(dispersion = "character",
                        fit_full = "list"))

################################################################################

setMethod("show", "dmSQTLfit", function(object){
  
  callNextMethod(object)
  
})


################################################################################

#' @rdname dmFit
#' @export
setMethod("dmFit", "dmSQTLdispersion", function(x, dispersion = "genewise_dispersion", prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  
  fit_full <- dmSQTL_fitOneModel(counts = x@counts, genotypes = x@genotypes, dispersion = slot(x, dispersion), model = "full", prop_mode = prop_mode, prop_tol = prop_tol, verbose = TRUE, BPPARAM = BPPARAM)
  
  
  return(new("dmSQTLfit", dispersion = dispersion, fit_full = fit_full, mean_expression = x@mean_expression, common_dispersion = x@common_dispersion, genewise_dispersion = x@genewise_dispersion, counts = x@counts, genotypes = x@genotypes, blocks = x@blocks, samples = x@samples))
  
  
})



################################################################################

#' @param snp_id Vector of SNP IDs to be plotted. \code{snp_id} must match \code{gene_id}.
#' @rdname plotFit
#' @export
setMethod("plotFit", "dmSQTLfit", function(x, gene_id, snp_id, plot_type = "boxplot1", order = TRUE, plot_full = TRUE, plot_main = TRUE, out_dir = NULL){
  
  stopifnot(plot_type %in% c("barplot", "boxplot1", "boxplot2", "lineplot", "ribbonplot"))
  
  
  dmSQTL_plotFit(gene_id = gene_id, snp_id = snp_id, counts = x@counts, genotypes = x@genotypes, blocks = x@blocks, samples = x@samples, dispersion = slot(x, x@dispersion), fit_full = x@fit_full, fit_null = NULL, table = NULL, plot_type = plot_type, order = order, plot_full = plot_full, plot_null = FALSE, plot_main = plot_main, out_dir = out_dir)
  
})














































