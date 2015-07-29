setClass("dmSQTLdata", representation(counts = "list", genotypes = "list", samples = "data.frame"))

dmSQTLdata <- function(counts, gene_id_counts, feature_id_counts, genotypes, gene_id_genotypes, snp_id_genotypes, sample_id){
  
  stopifnot( class( counts ) %in% c("matrix"))
  counts <- ceiling(as.matrix(counts))
  
  stopifnot( class( gene_id_counts ) %in% c("character", "factor"))
  stopifnot( class( feature_id_counts ) %in% c("character", "factor"))
  stopifnot( length(gene_id_counts) == length( feature_id_counts ) )

  stopifnot( class( genotypes ) %in% c("matrix"))
  
  stopifnot( class( gene_id_genotypes ) %in% c("character", "factor"))
  stopifnot( class( snp_id_genotypes ) %in% c("character", "factor"))
  stopifnot( length(gene_id_genotypes) == length( snp_id_genotypes ) )
  
  stopifnot( class( sample_id ) %in% c("character", "factor"))


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



dmSQTLdataFromRanges <- function(counts, gene_id_counts, feature_id_counts, gene_ranges, genotypes, snp_id_genotypes, snp_ranges, sample_id, window = 5e3){

  rownames(genotypes) <- snp_id_genotypes
  gene_ranges <- GenomicRanges::resize(gene_ranges, GenomicRanges::width(gene_ranges) + 2 * window, fix = "center")

  ## Match genes and SNPs
  variantMatch <- GenomicRanges::findOverlaps(gene_ranges, snp_ranges, select = "all")

  q <- queryHits(variantMatch)
  s <- subjectHits(variantMatch)

  genotypes_match <- genotypes[s, ]
  snp_id_genotypes_match <- snp_id_genotypes[s]
  gene_id_genotypes_match <- names(gene_ranges)[q]


  data <- dmSQTLdata(counts = counts, gene_id_counts = gene_id_counts, feature_id_counts = feature_id_counts, genotypes = genotypes_match, gene_id_genotypes = gene_id_genotypes_match, snp_id_genotypes = snp_id_genotypes_match, sample_id = sample_id)

  return(data)

}


