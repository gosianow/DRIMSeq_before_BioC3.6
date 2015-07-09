setClass("dmDSdata", representation(counts = "list", samples = "data.frame"))

dmDSdata <- function(counts, gene_id_counts, feature_id_counts, sample_id, group){
	
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



setClass("dmSQTLdata", representation(counts = "list", genotypes = "list", samples = "data.frame"))

dmSQTLdata <- function(counts, gene_id_counts, feature_id_counts, genotypes, gene_id_genotypes, snp_id_genotypes, sample_id){
	
	stopifnot( class( counts ) %in% c("matrix", "data.frame"))
	counts <- ceiling(as.matrix(counts))
	
	stopifnot( class( gene_id_counts ) %in% c("character", "factor"))
	stopifnot( class( feature_id ) %in% c("character", "factor"))
	stopifnot( class( sample_id ) %in% c("character", "factor"))
	stopifnot( length(gene_id_counts) == length( feature_id ) )
	
	genes2keep <- gene_id_counts %in% gene_id_genotypes
	counts <- counts[genes2keep, , drop = FALSE]
	gene_id_counts <- gene_id_counts[genes2keep]
	feature_id_counts <- feature_id_counts[genes2keep]

	genes2keep <- gene_id_genotypes %in% gene_id_counts
	genotypes <- genotypes[genes2keep, , drop = FALSE]
	gene_id_genotypes <- gene_id_genotypes[genes2keep]
	snp_id_genotypes <- snp_id_genotypes[genes2keep]


	colnames(counts) <- sample_id
	rownames(counts) <- feature_id_counts
	counts <- split.data.frame(counts, factor(gene_id_counts, levels = unique(gene_id_counts)))
	
	colnames(genotypes) <- sample_id
	rownames(genotypes) <- snp_id_genotypes
	genotypes <- split.data.frame(genotypes, factor(gene_id_genotypes, levels = unique(gene_id_counts)))

	samples <- data.frame(sample_id = sample_id)
	
	data <- new("dmSQTLdata", counts = counts, genotypes = genotypes, samples = samples)
	
	return(data)
	
}




















