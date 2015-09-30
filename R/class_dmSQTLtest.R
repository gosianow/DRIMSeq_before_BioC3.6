#' @include class_dmSQTLfit.R class_dmDStest.R
NULL

################################################################################

#' Object that extends \code{dmSQTLfit} by adding the test results.
#' 
#' @slot fit_null list
#' @slot results data.frame with \code{gene_id} - gene IDs, \code{snp_id} - SNP IDs, \code{lr} - likelihood ratio statistics, \code{df} - degrees of freedom, \code{pvalue} - p-values and \code{adj_pvalue} - Benjamini & Hochberg adjusted p-values.
setClass("dmSQTLtest", 
         contains = "dmSQTLfit",
         representation(fit_null = "list",
          results = "data.frame"))

##############################################################

#' @rdname dmSQTLtest-class
#' @export
setMethod("results", "dmSQTLtest", function(x) x@results )



################################################################################

setMethod("show", "dmSQTLtest", function(object){
  
  callNextMethod(object)
  
  cat("* data accessors: results()\n")
  
})


################################################################################
# prop_mode = "constrOptimG"; prop_tol = 1e-12; verbose = FALSE; BPPARAM = BiocParallel::MulticoreParam(workers = 10)

#' @rdname dmTest
#' @export
setMethod("dmTest", "dmSQTLfit", function(x, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  fit_null <- dmSQTL_fitOneModel(counts = x@counts, genotypes = x@genotypes, dispersion = slot(x, x@dispersion), model = "null", prop_mode = prop_mode, prop_tol = prop_tol, verbose = TRUE, BPPARAM = BPPARAM)
  
  results <- dmSQTL_test(fit_full = x@fit_full, fit_null = fit_null, BPPARAM = BPPARAM)
  
  colnames(results)[colnames(results) == "snp_id"] <- "block_id" 
  
  
  results_spl <- split(results, factor(results$gene_id, levels = names(x@blocks)))
  
  inds <- 1:length(results_spl)
  
  results_new <- lapply(inds, function(i){
    # i = 1
    
    res <- results_spl[[i]]
    blo <- x@blocks[[i]]
    
    matching <- match(blo[, "block_id"], res[, "block_id"])
    snp_id <- blo[, "snp_id"]
    
    res_new <- cbind(res[matching, c("gene_id", "block_id")], snp_id, res[matching, c("lr", "df", "pvalue", "adj_pvalue")])
    
    return(res_new)

    })
  
  results_new <- do.call(rbind, results_new)
  
  
  return(new("dmSQTLtest", fit_null = fit_null, results = results_new, dispersion = x@dispersion, fit_full = x@fit_full, mean_expression = x@mean_expression, common_dispersion = x@common_dispersion, genewise_dispersion = x@genewise_dispersion, counts = x@counts, genotypes = x@genotypes, blocks = x@blocks, samples = x@samples))
  
  
})


################################################################################

#' @rdname plotTest
#' @export
setMethod("plotTest", "dmSQTLtest", function(x, out_dir = NULL){
  
  dm_plotTable(pvalues = unique(x@results[, c("gene_id", "block_id", "pvalue")])[, "pvalue"], out_dir = out_dir)
  
})



################################################################################

#' @rdname plotFit
#' @export
setMethod("plotFit", "dmSQTLtest", function(x, gene_id, snp_id, plot_type = "boxplot1", order = TRUE, plot_full = TRUE, plot_null = TRUE, plot_main = TRUE, out_dir = NULL){
  
  stopifnot(plot_type %in% c("barplot", "boxplot1", "boxplot2", "lineplot", "ribbonplot"))
  
  dmSQTL_plotFit(gene_id = gene_id, snp_id = snp_id, counts = x@counts, genotypes = x@genotypes, blocks = x@blocks, samples = x@samples, dispersion = slot(x, x@dispersion), fit_full = x@fit_full, fit_null = x@fit_null, table = x@results, plot_type = plot_type, order = order, plot_full = plot_full, plot_null = plot_null, plot_main = plot_main, out_dir = out_dir)
  
  
})













































