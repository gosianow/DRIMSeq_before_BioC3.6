setClass("dmDSfit", 
  contains = "dmDSdispersion",
  representation(fit_full = "MatrixList",
    fit_null = "MatrixList"))


setMethod("show", "dmDSfit", function(object){
  
  callNextMethod(object)
  
  cat("Slot \"fit_full\":\n")
  print(object@fit_full)
  
  cat("Slot \"fit_null\":\n")
  print(object@fit_null)
  
  })



setGeneric("dmDSfit", function(x, ...) standardGeneric("dmDSfit"))


setMethod("dmDSfit", "dmDSdispersion", function(x, dispersion = "tagwise_dispersion", prop_mode = c("constrOptim", "constrOptimG", "FisherScoring")[2], prop_tol = 1e-12, verbose = FALSE, BPPARAM = MulticoreParam(workers=1)){
   
   
  
   
  fit_full <- dmDS_fitOneModel(counts = x@counts, samples = x@samples, dispersion = slot(x, dispersion), model = "full", prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)  
   
  fit_null <- dmDS_fitOneModel(counts = x@counts, samples = x@samples, dispersion = slot(x, dispersion), model = "null", prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
  
  
  
  return(new("dmDSfit", mean_expression = mean_expression, common_dispersion = common_dispersion, tagwise_dispersion = tagwise_dispersion, counts = x@counts, samples = x@samples))
  
  
  })




















































