setClass("dmDSdata", representation(counts = "list", samples = "data.frame"))

dmDSdata <- function(counts, gene_id, feature_id, sample_id, group){
	
	stopifnot( class( counts ) %in% c("matrix", "data.frame"))
	counts <- ceiling(as.matrix(counts))
	
  stopifnot( class( gene_id ) %in% c("character", "factor"))
  stopifnot( class( feature_id ) %in% c("character", "factor"))
  stopifnot( class( sample_id ) %in% c("character", "factor"))
  stopifnot( class( group ) %in% c("character", "factor"))
	stopifnot( length(gene_id) == length( feature_id ) )
	stopifnot( length(sample_id) == length( group ) )
	
	
	colnames(counts) <- sample_id
	rownames(counts) <- feature_id	
	counts <- split.data.frame(counts, factor(gene_id, levels = unique(gene_id)))
	
	samples <- data.frame(sample_id = sample_id, group = group)
	
	data <- new("dmDSdata", counts = counts, samples = samples)
	
	return(data)
	
}







