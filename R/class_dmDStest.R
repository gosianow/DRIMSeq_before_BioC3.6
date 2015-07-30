setClass("dmDStest", 
  contains = "dmDSfit",
  representation(table = "DataFrame"))


setMethod("show", "dmDStest", function(object){
  
  callNextMethod(object)
  
  cat("Slot \"table\":\n")
  print(object@table)
  
  
  })



setGeneric("dmDStest", function(x, ...) standardGeneric("dmDStest"))


setMethod("dmDStest", "dmDSfit", function(x){
   
   
  table <- dmDS_test(stats_full = x@statistics_full, stats_null = x@statistics_null)
  
  
  return(new("dmDStest", table = table, dispersion = x@dispersion, proportions_full = x@proportions_full, statistics_full = x@statistics_full, proportions_null = x@proportions_null, statistics_null = x@statistics_null, mean_expression = x@mean_expression, common_dispersion = x@common_dispersion, tagwise_dispersion = x@tagwise_dispersion, counts = x@counts, samples = x@samples))
  
      
  })



setGeneric("dmDSplotTest", function(x, ...) standardGeneric("dmDSplotTest"))

setMethod("dmDSplotTest", "dmDStest", function(x, out_dir = NULL){
  
  dm_plotTable(x@table, out_dir = out_dir)
  
  })



setMethod("dmDSplotFit", "dmDStest", function(x, gene_id, plot_type = c("barplot", "boxplot1", "boxplot2", "lineplot", "ribbonplot")[1], order = TRUE, plot_full = TRUE, plot_nunll = TRUE, out_dir = NULL){
  
  
  dmDS_plotFit(gene_id = gene_id, counts = x@counts, samples = x@samples, dispersion = slot(x, x@dispersion), proportions_full = x@proportions_full, proportions_null = x@proportions_null, table = x@table, plot_type = plot_type, order = order, plot_full = plot_full, plot_null = plot_null, out_dir = out_dir)
  
  
  })













































