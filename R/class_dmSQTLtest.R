setClass("dmSQTLtest", 
  contains = "dmSQTLfit",
  representation(table = "DataFrame"))


setMethod("show", "dmSQTLtest", function(object){
  
  callNextMethod(object)
  
  cat("\nSlot \"table\":\n")
  print(object@table)
  
  
  })



setGeneric("dmSQTLtest", function(x, ...) standardGeneric("dmSQTLtest"))


setMethod("dmSQTLtest", "dmSQTLfit", function(x, BPPARAM = MulticoreParam(workers=1)){
  
  # fit_full = x@fit_full; fit_null = x@fit_null
  
  table <- dmSQTL_test(fit_full = x@fit_full, fit_null = x@fit_null, BPPARAM = BPPARAM)
  
  
  return(new("dmSQTLtest", table = table, dispersion = x@dispersion, fit_full = x@fit_full, fit_null = x@fit_null, mean_expression = x@mean_expression, common_dispersion = x@common_dispersion, tagwise_dispersion = x@tagwise_dispersion, counts = x@counts, genotypes = x@genotypes, samples = x@samples))
  
  
  })



setGeneric("dmSQTLplotTest", function(x, ...) standardGeneric("dmSQTLplotTest"))

setMethod("dmSQTLplotTest", "dmSQTLtest", function(x, out_dir = NULL){
  
  dm_plotTable(x@table, out_dir = out_dir)
  
  })



setMethod("dmSQTLplotFit", "dmSQTLtest", function(x, gene_id, snp_id, plot_type = c("barplot", "boxplot1", "boxplot2", "lineplot", "ribbonplot")[1], order = TRUE, plot_full = TRUE, plot_nunll = TRUE, out_dir = NULL){
  
  
  dmSQTL_plotFit(gene_id = gene_id, snp_id = snp_id, counts = x@counts, genotypes = x@genotypes, samples = x@samples, dispersion = slot(x, x@dispersion), fit_full = x@fit_full, fit_null = x@fit_null, table = x@table, plot_type = plot_type, order = order, plot_full = plot_full, plot_null = plot_null, out_dir = out_dir)
 
 
 })













































