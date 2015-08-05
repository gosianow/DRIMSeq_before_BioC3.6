#' @include class_dmFit.R
NULL

setClass("dmDSfit", 
  contains = "dmDSdispersion",
  representation(dispersion = "character",
    fit_full = "dmFit",
    fit_null = "dmFit"))


setMethod("show", "dmDSfit", function(object){
  
  callNextMethod(object)
  
  cat("\nSlot \"dispersion\":\n")
  print(object@dispersion)
  
  cat("\nSlot \"fit_full\":\n")
  print(object@fit_full)
  
  cat("\nSlot \"fit_null\":\n")
  print(object@fit_null)
  
  
  })



setGeneric("dmDSfit", function(x, ...) standardGeneric("dmDSfit"))


# counts = x@counts; samples = x@samples; dispersion = "tagwise_dispersion"; prop_mode = c("constrOptim", "constrOptimG", "FisherScoring")[2]; prop_tol = 1e-12; verbose = FALSE; BPPARAM = MulticoreParam(workers=15)

setMethod("dmDSfit", "dmDSdispersion", function(x, dispersion = "tagwise_dispersion", prop_mode = c("constrOptim", "constrOptimG", "FisherScoring")[2], prop_tol = 1e-12, verbose = FALSE, BPPARAM = MulticoreParam(workers=1)){
   
   
  fit_full <- dmDS_fitOneModel(counts = x@counts, samples = x@samples, dispersion = slot(x, dispersion), model = "full", prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)

   
  fit_null <- dmDS_fitOneModel(counts = x@counts, samples = x@samples, dispersion = slot(x, dispersion), model = "null", prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
  

  return(new("dmDSfit", dispersion = dispersion, fit_full = fit_full, fit_null = fit_null, mean_expression = x@mean_expression, common_dispersion = x@common_dispersion, tagwise_dispersion = x@tagwise_dispersion, counts = x@counts, samples = x@samples))
  
  
  })


setGeneric("dmDSplotFit", function(x, ...) standardGeneric("dmDSplotFit"))

setMethod("dmDSplotFit", "dmDSfit", function(x, gene_id, plot_type = c("barplot", "boxplot1", "boxplot2", "lineplot", "ribbonplot")[1], order = TRUE, plot_full = TRUE, plot_nunll = TRUE, out_dir = NULL){
  
  
  dmDS_plotFit(gene_id = gene_id, counts = x@counts, samples = x@samples, dispersion = slot(x, x@dispersion), proportions_full = x@fit_full@proportions, proportions_null = x@fit_null@proportions, table = NULL, plot_type = plot_type, order = order, plot_full = plot_full, plot_null = plot_null, out_dir = out_dir)
  
  
  })

















































