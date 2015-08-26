#' @include class_dmDSfit.R
NULL

################################################################################

#' Object that extends \code{dmDSfit} by adding the test results.
#' 
#' @slot fit_null list of null models defined by contrasts.
#' @slot results data.frame with \code{gene_id} - gene IDs, \code{lr} - likelihood ratio statistics, \code{df} - degrees of freedom, \code{pvalue} - p-values and \code{adj_pvalue} - Benjamini & Hochberg adjusted p-values.
setClass("dmDSLRT", 
         contains = "dmDSfit",
         representation(pairwise_comparison = "matrix",
          fit_null = "list",
          results = "data.frame"))

################################################################################

setMethod("show", "dmDSLRT", function(object){
  
  callNextMethod(object)
  
  cat("\nSlot \"pairwise_comparison\":\n")
  show_matrix(object@pairwise_comparison, nhead = 5, ntail = 5)
  
  cat("\nSlot \"fit_null\":\n")
  show_MatrixList_list(fit_null)
  
  cat("\nSlot \"results\":\n")
  show_matrix(object@results, nhead = 5, ntail = 5)
  
  
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
#' @inheritParams dmFit
#' @param pairwise_comparison matrix
#' @return This function returns a \code{\linkS4class{dmDSLRT}} or \code{\linkS4class{dmSQTLLRT}} object with an additional slot \code{table} which is sorted by significance and contains  \code{gene_id} - gene IDs, \code{lr} - likelihood ratio statistics, \code{df} - degrees of freedom, \code{pvalue} - p-values and \code{adj_pvalue} - Benjamini & Hochberg adjusted p-values.
#' @export
setMethod("dmLRT", "dmDSfit", function(x, pairwise_comparison = matrix(nrow = 0, ncol = 2), prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  np <- nrow(pairwise_comparison)
  
  if(np == 0){
    message("Running comparison between all groups..")
    fit_null <- list()
    
    fit_null[[1]] <- dmDS_fitOneModel(counts = x@counts, samples = x@samples, dispersion = slot(x, x@dispersion), model = "null", prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
    
    results <- dmDS_test(stats_full = x@fit_full@metadata, stats_null = fit_null[[1]]@metadata)
    
    
  }else{
    
    fit_null <- vector("list", np)
    names(fit_null) <- rownames(pairwise_comparison)
    tables <- vector("list", np)
    
    if(is.null(rownames(pairwise_comparison)))
    suffix <- 1:np
    else
    suffix <- rownames(pairwise_comparison)
    
    
    for(i in 1:np){
      # i = 1
      message("Running comparison between group ", pairwise_comparison[i,1], " and group ", pairwise_comparison[i,2], "..")
      
      if(mode(pairwise_comparison) == "numeric")
      samps <- as.numeric(x@samples$group) %in% pairwise_comparison[i, ]
      if(mode(pairwise_comparison) == "character")
      samps <- x@samples$group %in% pairwise_comparison[i, ]
      
      samples = x@samples[samps, ]
      samples$sample_id <- factor(samples$sample_id)
      samples$group <- factor(samples$group)
      
      fit_null[[i]] <- dmDS_fitOneModel(counts = x@counts[, samps], samples = samples, dispersion = slot(x, x@dispersion), model = "null", prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
      
      tables[[i]] <- dmDS_test(stats_full = x@fit_full@metadata[, pairwise_comparison[i, ]], stats_null = fit_null[[i]]@metadata)
      
      if(np > 1)
      names(tables[[i]])[-1] <- paste0(names(tables[[i]])[-1], "_", suffix[i])
      
      
    }
    
    if(np > 1)
    results <- Reduce(function(...) merge(..., by = "gene_id", all = TRUE, sort = FALSE), tables)
    else
    results <- tables[[1]]
    
  }
  
  return(new("dmDSLRT", pairwise_comparison = pairwise_comparison, fit_null = fit_null, results = results, dispersion = x@dispersion, fit_full = x@fit_full,  mean_expression = x@mean_expression, common_dispersion = x@common_dispersion, genewise_dispersion = x@genewise_dispersion, counts = x@counts, samples = x@samples))
  
  
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
  
  dm_plotTable(x@results, out_dir = out_dir)
  
})

################################################################################

#' @rdname plotFit
#' @export
setMethod("plotFit", "dmDSLRT", function(x, gene_id, plot_type = c("barplot", "boxplot1", "boxplot2", "lineplot", "ribbonplot")[1], order = TRUE, plot_full = TRUE, plot_null = TRUE, out_dir = NULL){
  
  stopifnot(plot_type %in% c("barplot", "boxplot1", "boxplot2", "lineplot", "ribbonplot"))
  
  dmDS_plotFit(gene_id = gene_id, counts = x@counts, samples = x@samples, dispersion = slot(x, x@dispersion), proportions_full = x@fit_full, proportions_null = x@fit_full, table = x@table, plot_type = plot_type, order = order, plot_full = plot_full, plot_null = plot_null, out_dir = out_dir)
  
  
})













































