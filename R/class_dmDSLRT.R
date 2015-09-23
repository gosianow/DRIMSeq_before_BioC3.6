#' @include class_dmDSfit.R
NULL

################################################################################

#' Object that extends \code{dmDSfit} by adding the test results.
#' 
#' @slot fit_null list of null models defined by contrasts.
#' @slot results data.frame with \code{gene_id} - gene IDs, \code{lr} - likelihood ratio statistics, \code{df} - degrees of freedom, \code{pvalue} - p-values and \code{adj_pvalue} - Benjamini & Hochberg adjusted p-values.
setClass("dmDSLRT", 
         contains = "dmDSfit",
         representation(compared_groups = "list",
          fit_null = "list",
          results = "data.frame"))


##############################################################

#' @rdname dmDSLRT-class
#' @export
setMethod("proportions", "dmDSLRT", function(x){
  
  nc <- length(x@compared_groups)
  
  if(nc == 1){
    
    prop_null <- x@fit_null@unlistData
    
  }else{
    
    prop_null <- do.call(cbind, lapply(x@fit_null, function(f) f@unlistData))
    
    if(is.null(names(x@compared_groups)))
    suffix <- 1:nc
    else
    suffix <- names(x@compared_groups)
    
    colnames(prop_null) <- paste0("null_", suffix)
    
  } 
  
  data.frame(gene_id = rep(names(x@counts), width(x@counts)), feature_id = rownames(x@counts@unlistData), x@fit_full@unlistData, prop_null, stringsAsFactors = FALSE, row.names = NULL)
  
})


#' @rdname dmDSLRT-class
#' @export
setMethod("statistics", "dmDSLRT", function(x){
  
  nc <- length(x@compared_groups)
  
  if(nc == 1){
    
    stats_null <- x@fit_null@metadata
    
  }else{
    
    stats_null <- do.call(cbind, lapply(x@fit_null, function(f) f@metadata))
    
    if(is.null(names(x@compared_groups)))
    suffix <- 1:nc
    else
    suffix <- names(x@compared_groups)
    
    colnames(stats_null) <- paste0(colnames(stats_null), "_null_", rep(suffix, each = 2))
    
  } 
  
  
  df <- data.frame(gene_id = names(x@counts), x@fit_full@metadata, stats_null, stringsAsFactors = FALSE, row.names = NULL)
  colnames(df)[2:(ncol(x@fit_full@metadata)+1)] <- paste0("lik_", colnames(x@fit_full@metadata))
  
  return(df)
  
})



#' @rdname dmDSLRT-class
#' @export
setGeneric("results", function(x, ...) standardGeneric("results"))

#' @rdname dmDSLRT-class
#' @export
setMethod("results", "dmDSLRT", function(x) x@results )



################################################################################

setMethod("show", "dmDSLRT", function(object){
  
  callNextMethod(object)
  
  cat("  results()\n")
  
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
#' @param compared_groups vector or a list of vectors.
#' @return This function returns a \code{\linkS4class{dmDSLRT}} or \code{\linkS4class{dmSQTLLRT}} object with an additional slot \code{table} which is sorted by significance and contains  \code{gene_id} - gene IDs, \code{lr} - likelihood ratio statistics, \code{df} - degrees of freedom, \code{pvalue} - p-values and \code{adj_pvalue} - Benjamini & Hochberg adjusted p-values.
#' @examples 
#' d <- dataDS_dmDSdispersion
#' 
#' # If possible, increase the number of workers
#' d <- dmFit(d, BPPARAM = BiocParallel::MulticoreParam(workers = 1))
#' 
#' d <- dmLRT(d)
#' 
#' results <- results(d)
#' 
#' @export
setMethod("dmLRT", "dmDSfit", function(x, compared_groups = 1:nlevels(samples(x)$group), prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  if(!is.list(compared_groups))
  compared_groups <- list(compared_groups)
  
  nc <- length(compared_groups)
  
    fit_null <- vector("list", nc)
    names(fit_null) <- names(compared_groups)
    tables <- vector("list", nc)
    
    suffix <- NULL
    
    if(nc > 1){
      if(is.null(names(compared_groups)))
      suffix <- 1:nc
      else
      suffix <- names(compared_groups)
    }
    
    
    for(i in 1:nc){
      # i = 1
      
      if(mode(compared_groups[[i]]) == "numeric")
      samps <- as.numeric(x@samples$group) %in% compared_groups[[i]]
      if(mode(compared_groups[[i]]) == "character")
      samps <- x@samples$group %in% compared_groups[[i]]
      
      samples = x@samples[samps, ]
      samples$sample_id <- factor(samples$sample_id)
      samples$group <- factor(samples$group)
      
      message("Running comparison ", suffix[i], " between groups: ", paste0(levels(samples$group), collapse = ", "))
      
      
      fit_null[[i]] <- dmDS_fitOneModel(counts = x@counts[, samps], samples = samples, dispersion = slot(x, x@dispersion), model = "null", prop_mode = prop_mode, prop_tol = prop_tol, verbose = TRUE, BPPARAM = BPPARAM)
      
      tables[[i]] <- dmDS_test(stats_full = x@fit_full@metadata[, compared_groups[[i]]], stats_null = fit_null[[i]]@metadata)
      
      if(nc > 1)
      names(tables[[i]])[-1] <- paste0(names(tables[[i]])[-1], "_", suffix[i])
      
      
    }
    
    if(nc > 1)
    results <- Reduce(function(...) merge(..., by = "gene_id", all = TRUE, sort = FALSE), tables)
    else
    results <- tables[[1]]
    
  
  return(new("dmDSLRT", compared_groups = compared_groups, fit_null = fit_null, results = results, dispersion = x@dispersion, fit_full = x@fit_full,  mean_expression = x@mean_expression, common_dispersion = x@common_dispersion, genewise_dispersion = x@genewise_dispersion, counts = x@counts, samples = x@samples))
  
  
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
#' @examples 
#' d <- dataDS_dmDSdispersion
#' 
#' # If possible, increase the number of workers
#' d <- dmFit(d, BPPARAM = BiocParallel::MulticoreParam(workers = 1))
#' 
#' d <- dmLRT(d)
#' 
#' results <- results(d)
#' 
#' plotLRT(d)
#' plot(d)
#' 
#' @export
setMethod("plotLRT", "dmDSLRT", function(x, out_dir = NULL){
  
  col_pv <- colnames(x@results)[grepl("pvalue", colnames(x@results)) & !grepl("adj_pvalue", colnames(x@results))]
  
  for(i in 1:length(col_pv))
  dm_plotTable(pvalues = x@results[, col_pv[i]], name = col_pv[i], out_dir = out_dir)
  
})


##############################################################

#' @rdname dmDSLRT-class
#' @export
setMethod("plot", "dmDSLRT", function(x, out_dir = NULL){
  
  plotLRT(x, out_dir = out_dir)
  
})



################################################################################

#' @rdname plotFit
#' @param compared_groups numeric or character indicating the comparison that should be plotted.
#' @export
setMethod("plotFit", "dmDSLRT", function(x, gene_id, plot_type = "barplot", order = TRUE, plot_full = TRUE, plot_null = TRUE, compared_groups = 1, plot_main = TRUE, out_dir = NULL){
  
  stopifnot(plot_type %in% c("barplot", "boxplot1", "boxplot2", "lineplot", "ribbonplot"))
  
  i <- compared_groups
  compared_groups <- x@compared_groups
  
  
  if(mode(compared_groups[[i]]) == "numeric")
  samps <- as.numeric(x@samples$group) %in% compared_groups[[i]]
  if(mode(compared_groups[[i]]) == "character")
  samps <- x@samples$group %in% compared_groups[[i]]
  
  samples = x@samples[samps, ]
  samples$sample_id <- factor(samples$sample_id)
  samples$group <- factor(samples$group)
  
  
  if(length(compared_groups) > 1){
    
    if(is.null(names(compared_groups))){
      suffix <- i
    }else{
      if(mode(i) == "numeric")
      suffix <- names(compared_groups)[i]
      else
      suffix <- i
    }
    
    which_cols <- grepl(paste0("_", suffix), colnames(x@results)) | grepl("gene_id", colnames(x@results)) 
    results <- x@results[, which_cols]
    colnames(results) <- gsub(paste0("_", suffix), "", colnames(results))
    
  }else{
    results <- x@results
  }

  
  dmDS_plotFit(gene_id = gene_id, counts = x@counts[, samps], samples = samples, dispersion = slot(x, x@dispersion), proportions_full = x@fit_full[, compared_groups[[i]]], proportions_null = x@fit_null[[i]], table = results, plot_type = plot_type, order = order, plot_full = plot_full, plot_null = plot_null, plot_main = plot_main, out_dir = out_dir)
  
  
})













































