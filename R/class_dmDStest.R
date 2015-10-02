#' @include class_dmDSfit.R
NULL

################################################################################

#' dmDStest object
#' 
#' dmDStest extends the \code{\linkS4class{dmDSfit}} class by adding the null model Dirichlet-multinomial feature proportion estimates and the results of testing for differential splicing. Result of \code{\link{dmTest}}.
#' 
#' @details 
#' 
#' \itemize{
#'  \item \code{proportions(x)}: Get a data.frame with estimated feature ratios for full model and null models specified in \code{\link{dmTest}} with \code{compared_groups} parameter.
#'  \item \code{statistics(x)}: Get a data.frame with full and null log-likelihoods and degrees of freedom.
#'  \item \code{results(x)}: Get a data.frame with results.
#' }
#' 
#' @param x dmDStest object.
#' @param ... Other parameters that can be defined by methods using this generic.
#' 
#' @slot compared_groups List of vectors specifying which group comparisons are performed.
#' @slot fit_null List of \code{\linkS4class{MatrixList}}. Each of them contains null proportions, likelihoods and degrees of freedom for a comparison specified with \code{compared_groups}.
#' @slot results data.frame with \code{gene_id} - gene IDs, \code{lr} - likelihood ratio statistics, \code{df} - degrees of freedom, \code{pvalue} - p-values and \code{adj_pvalue} - Benjamini & Hochberg adjusted p-values.
#' @author Malgorzata Nowicka
#' @seealso \code{\link{plotTest}}, \code{\linkS4class{dmDSdata}}, \code{\linkS4class{dmDSdispersion}}, \code{\linkS4class{dmDSfit}}
setClass("dmDStest", 
         contains = "dmDSfit",
         representation(compared_groups = "list",
          fit_null = "list",
          results = "data.frame"))


##############################################################

#' @rdname dmDStest-class
#' @export
setMethod("proportions", "dmDStest", function(x){
  
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


#' @rdname dmDStest-class
#' @export
setMethod("statistics", "dmDStest", function(x){
  
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



#' @rdname dmDStest-class
#' @export
setGeneric("results", function(x, ...) standardGeneric("results"))

#' @rdname dmDStest-class
#' @export
setMethod("results", "dmDStest", function(x) x@results )



################################################################################

setMethod("show", "dmDStest", function(object){
  
  callNextMethod(object)
  
  cat("  results()\n")
  
})

################################################################################

#' Likelihood ratio test
#' 
#' First, estimate the null Dirichlet-multinomial model proportions, i.e., feature ratios are estimated based on pooled data. Use the likelihood ratio statistic to test for the difference between feature proportions in different groups to identify the differentially spliced genes (differential splicing analysis) or the sQTLs (sQTL analysis).
#' 
#' @param x \code{\linkS4class{dmDSfit}} or \code{\linkS4class{dmSQTLfit}} object.
#' @param ... Other parameters that can be defined by methods using this generic.
#' @export
setGeneric("dmTest", function(x, ...) standardGeneric("dmTest"))


################################################################################


#' @inheritParams dmFit
#' @param compared_groups Vector or a list of vectors that defines which experimental conditions should be tested for differential splicing. By default, we test for a difference between any of the groups specified in \code{samples(x)$group}. Values in this vectors should be levels or numbers of the levels of \code{samples(x)$group}.
#' 
#' @return Returns a \code{\linkS4class{dmDStest}} or \code{\linkS4class{dmSQTLtest}} object.
#' @examples 
#' ### Differential splicing analysis
#' 
#' d <- dataDS_dmDSdispersion
#' 
#' # If possible, increase the number of workers
#' d <- dmFit(d, BPPARAM = BiocParallel::MulticoreParam(workers = 1))
#' 
#' d <- dmTest(d, BPPARAM = BiocParallel::MulticoreParam(workers = 1))
#' 
#' plotTest(d)
#' 
#' res <- results(d)
#' res <- res[order(res$pvalue, decreasing = FALSE), ]
#' 
#' plotFit(d, gene_id = res$gene_id[1])
#' 
#' @author Malgorzata Nowicka
#' @seealso \code{\link{plotTest}}, \code{\link{dmDispersion}}, \code{\link{dmFit}}
#' @rdname dmTest
#' @export
setMethod("dmTest", "dmDSfit", function(x, compared_groups = 1:nlevels(samples(x)$group), prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
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
    
  
  return(new("dmDStest", compared_groups = compared_groups, fit_null = fit_null, results = results, dispersion = x@dispersion, fit_full = x@fit_full,  mean_expression = x@mean_expression, common_dispersion = x@common_dispersion, genewise_dispersion = x@genewise_dispersion, counts = x@counts, samples = x@samples))
  
  
})


################################################################################

#' Plot a histogram of p-values
#' 
#' @param x \code{\linkS4class{dmDStest}} or \code{\linkS4class{dmSQTLtest}} object.
#' @export
setGeneric("plotTest", function(x, ...) standardGeneric("plotTest"))



################################################################################

#' @inheritParams plotData
#' @examples
#' ### Differential splicing analysis
#' 
#' d <- dataDS_dmDStest
#' plotTest(d)
#' 
#' @author Malgorzata Nowicka
#' @seealso \code{\link{plotData}}, \code{\link{plotDispersion}}, \code{\link{plotFit}}
#' @rdname plotTest
#' @export
setMethod("plotTest", "dmDStest", function(x, out_dir = NULL){
  
  col_pv <- colnames(x@results)[grepl("pvalue", colnames(x@results)) & !grepl("adj_pvalue", colnames(x@results))]
  
  for(i in 1:length(col_pv))
  dm_plotTable(pvalues = x@results[, col_pv[i]], name = col_pv[i], out_dir = out_dir)
  
})



################################################################################

#' @param compared_groups Numeric or character indicating comparison that should be plotted.
#' @param plot_null Logical. Whether to plot the proportions estimated by the null model.
#' @rdname plotFit
#' @export
setMethod("plotFit", "dmDStest", function(x, gene_id, plot_type = "barplot", order = TRUE, plot_full = TRUE, plot_null = TRUE, compared_groups = 1, plot_main = TRUE, out_dir = NULL){
  
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













































