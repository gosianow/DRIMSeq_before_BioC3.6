#' @include class_dmFit.R
NULL


################################################################################
#' Object that extends \code{dmSQTLdispersion} by adding fitting.
#' 
#' @slot dispersion Character specifying which type of dispersion was used for fitting.
#' @slot fit_full A List of \code{\linkS4class{dmFit}} objects. Each of them contains the per gene fitting results that include the group proportions (groups are defined by genotypes) of feature expression, full model likelihoods and degrees of freedom.
#' @slot fit_null A List of \code{\linkS4class{dmFit}} objects. Each of them contsins the per gene fitting results that include the pooled proportions of feature expression, null model likelihoods and degrees of freedom.
#' @importClassesFrom S4Vectors List
setClass("dmSQTLfit", 
  contains = "dmSQTLdispersion",
  representation(dispersion = "character",
    fit_full = "List",
    fit_null = "List"))

################################################################################
setMethod("show", "dmSQTLfit", function(object){
  
  callNextMethod(object)
  
  cat("\nSlot \"dispersion\":\n")
  print(object@dispersion)
  
  cat("\nSlot \"fit_full\":\n")
  print(object@fit_full)
  
  cat("\nSlot \"fit_null\":\n")
  print(object@fit_null)
  
  
  })


################################################################################
#' Estimating the proportions and likelihoods of Dirichlet-multinomial full and null models.
#' 
#' @param x \code{\link{dmSQTLdispersion}} object or any that inherits from it i.e. \code{\link{dmSQTLfit}} or \code{\link{dmSQTLtest}}.
#' @param ... Parameters needed for the proportion estimation.
#' @export
setGeneric("dmSQTLfit", function(x, ...) standardGeneric("dmSQTLfit"))



################################################################################
#' @rdname dmSQTLfit
#' @inheritParams dmSQTLdispersion
#' @param dispersion Characted defining which dispersion should be used for fitting. Possible values \code{"tagwise_dispersion", "common_dispersion"}
#' @examples 
#' data <- dataSQTL_dmSQTLdispersion
#' data <- dmSQTLfit(data)
#' 
#' snp_id <- "snp_19_12796435"
#' gene_id <- "ENSG00000132004.7"
#' 
#' dmSQTLplotFit(data, gene_id, snp_id)
#' 
#' @export
setMethod("dmSQTLfit", "dmSQTLdispersion", function(x, dispersion = "tagwise_dispersion", prop_mode = c("constrOptim", "constrOptimG", "FisherScoring")[2], prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers=1)){
  
  # counts = x@counts; genotypes = x@genotypes; dispersion = slot(x, dispersion); model = "full"
   
  fit_full <- S4Vectors::List(dmSQTL_fitOneModel(counts = x@counts, genotypes = x@genotypes, dispersion = slot(x, dispersion), model = "full", prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM))
  
   
  fit_null <- S4Vectors::List(dmSQTL_fitOneModel(counts = x@counts, genotypes = x@genotypes, dispersion = slot(x, dispersion), model = "null", prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM))
  
 
  return(new("dmSQTLfit", dispersion = dispersion, fit_full = fit_full, fit_null = fit_null, mean_expression = x@mean_expression, common_dispersion = x@common_dispersion, tagwise_dispersion = x@tagwise_dispersion, counts = x@counts, genotypes = x@genotypes, samples = x@samples))
  
  
  })



################################################################################
#' Plot the estimated proportions.
#' 
#' @param x \code{\link{dmSQTLfit}} object or any that inherits from it i.e. \code{\link{dmSQTLtest}}.
#' @param ... Plotting parameters.
#' @export
setGeneric("dmSQTLplotFit", function(x, ...) standardGeneric("dmSQTLplotFit"))



################################################################################
#' @rdname dmSQTLplotFit
#' @inheritParams dmDSplotFit
#' @param snp_id Vector of SNP IDs to be plotted.
#' @export
setMethod("dmSQTLplotFit", "dmSQTLfit", function(x, gene_id, snp_id, plot_type = "boxplot1", order = TRUE, plot_full = TRUE, plot_null = TRUE, out_dir = NULL){
  
  # counts = x@counts; genotypes = x@genotypes; samples = x@samples; dispersion = slot(x, x@dispersion); fit_full = x@fit_full; fit_null = x@fit_null; table = NULL
  
  dmSQTL_plotFit(gene_id = gene_id, snp_id = snp_id, counts = x@counts, genotypes = x@genotypes, samples = x@samples, dispersion = slot(x, x@dispersion), fit_full = x@fit_full, fit_null = x@fit_null, table = NULL, plot_type = plot_type, order = order, plot_full = plot_full, plot_null = plot_null, out_dir = out_dir)
  
  
  })

















































