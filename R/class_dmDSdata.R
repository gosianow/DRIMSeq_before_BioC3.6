setClass("dmDSdata", representation(counts = "MatrixList", samples = "DataFrame"))

dmDSdata <- function(counts, gene_id_counts, feature_id_counts, sample_id, group){
  
  unlistData <- counts

  stopifnot(class(unlistData) == "matrix")
  unlistData <- ceiling(unlistData)
  
  stopifnot( class( gene_id_counts ) %in% c("character", "factor"))
  stopifnot( class( feature_id_counts ) %in% c("character", "factor"))
  stopifnot( class( sample_id ) %in% c("character", "factor"))
  stopifnot( class( group ) %in% c("character", "factor"))
  stopifnot( length(gene_id_counts) == length( feature_id_counts ) )
  stopifnot( length(sample_id) == length( group ) )
  stopifnot(nrow(unlistData) == length(feature_id_counts))
  stopifnot(ncol(unlistData) == length(sample_id))
  
  colnames(unlistData) <- sample_id
  rownames(unlistData) <- feature_id_counts
  
  counts <- new("MatrixList", unlistData = unlistData, partitioning = PartitioningByEnd(as.numeric(factor(gene_id_counts, levels = unique(gene_id_counts))), NG = length(unique(gene_id_counts)), names = unique(gene_id_counts)))
  
  samples <- DataFrame(sample_id = sample_id, group = group)
  
  data <- new("dmDSdata", counts = counts, samples = samples)
  
  return(data)
  
}

# setValidity("Greeting", function(object){
  
#   # has to return TRUE when valid object!
  
#   })


setMethod("show", "dmDSdata", function(object){
  
  cat("An object of class", class(object), "\n")
  
  cat("Slot \"counts\":\n")
  print(object@counts)
  
  cat("Slot \"samples\":\n")
  print(object@samples)
  
  })




setMethod("names", "dmDSdata", function(x) names(data@counts) )


setGeneric("dmDSfilter", function(x, ...) standardGeneric("dmDSfilter"))


setMethod("dmDSfilter", "dmDSdata", function(x, min_samps_gene_expr = 3, min_gene_expr = 1, min_samps_feature_prop = 3, min_feature_prop = 0.01, max_features = Inf){
  
  
  data_filtered <- dmDS_filter(counts = x@counts, samples = x@samples, min_samps_gene_expr = min_samps_gene_expr, min_gene_expr = min_gene_expr, min_samps_feature_prop = min_samps_feature_prop, min_feature_prop = min_feature_prop, max_features = max_features)
  
  return(data_filtered)
  
  
  })



setGeneric("dmDSplotData", function(x, ...) standardGeneric("dmDSplotData"))

setMethod("dmDSplotData", "dmDSdata", function(x, out_dir = NULL, info = NULL){
  
  dmDS_plotData(counts = x@counts, out_dir = out_dir, info = info)
  
  })


















































