setClass("dmSQTLfit", 
  contains = "dmSQTLdispersion",
  representation(dispersion = "character",
    fit_full = "List",
    fit_null = "List"))


setMethod("show", "dmSQTLfit", function(object){
  
  callNextMethod(object)
  
  cat("\nSlot \"dispersion\":\n")
  print(object@dispersion)
  
  cat("\nSlot \"fit_full\":\n")
  print(object@fit_full)
  
  cat("\nSlot \"fit_null\":\n")
  print(object@fit_null)
  
  
  })



setGeneric("dmSQTLfit", function(x, ...) standardGeneric("dmSQTLfit"))


setMethod("dmSQTLfit", "dmSQTLdispersion", function(x, dispersion = "tagwise_dispersion", prop_mode = c("constrOptim", "constrOptimG", "FisherScoring")[2], prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers=1)){
  
  # counts = x@counts; genotypes = x@genotypes; dispersion = slot(x, dispersion); model = "full"
   
  fit_full <- List(dmSQTL_fitOneModel(counts = x@counts, genotypes = x@genotypes, dispersion = slot(x, dispersion), model = "full", prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM))
  
   
  fit_null <- List(dmSQTL_fitOneModel(counts = x@counts, genotypes = x@genotypes, dispersion = slot(x, dispersion), model = "null", prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM))
  
 
  return(new("dmSQTLfit", dispersion = dispersion, fit_full = fit_full, fit_null = fit_null, mean_expression = x@mean_expression, common_dispersion = x@common_dispersion, tagwise_dispersion = x@tagwise_dispersion, counts = x@counts, genotypes = x@genotypes, samples = x@samples))
  
  
  })


setGeneric("dmSQTLplotFit", function(x, ...) standardGeneric("dmSQTLplotFit"))

setMethod("dmSQTLplotFit", "dmSQTLfit", function(x, gene_id, snp_id, plot_type = c("barplot", "boxplot1", "boxplot2", "lineplot", "ribbonplot")[1], order = TRUE, plot_full = TRUE, plot_nunll = TRUE, out_dir = NULL){
  
  # counts = x@counts; genotypes = x@genotypes; samples = x@samples; dispersion = slot(x, x@dispersion); fit_full = x@fit_full; fit_null = x@fit_null; table = NULL
  
  dmSQTL_plotFit(gene_id = gene_id, snp_id = snp_id, counts = x@counts, genotypes = x@genotypes, samples = x@samples, dispersion = slot(x, x@dispersion), fit_full = x@fit_full, fit_null = x@fit_null, table = NULL, plot_type = plot_type, order = order, plot_full = plot_full, plot_null = plot_null, out_dir = out_dir)
  
  
  })

















































