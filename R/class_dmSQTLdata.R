#' @importClassesFrom IRanges DataFrame
setClass("dmSQTLdata", 
  representation(counts = "MatrixList", 
    genotypes = "MatrixList", 
    samples = "DataFrame"))



dmSQTLdata <- function(counts, gene_id_counts, feature_id_counts, genotypes, gene_id_genotypes, snp_id_genotypes, sample_id){
  
  stopifnot( class( counts ) %in% c("matrix"))
  stopifnot( mode( counts ) %in% c("numeric"))
  counts <- ceiling(counts)
  
  stopifnot( class( gene_id_counts ) %in% c("character"))
  stopifnot( class( feature_id_counts ) %in% c("character"))
  stopifnot( length(gene_id_counts) == length( feature_id_counts ) )

  stopifnot( class( genotypes ) %in% c("matrix"))
  
  stopifnot( class( gene_id_genotypes ) %in% c("character"))
  stopifnot( class( snp_id_genotypes ) %in% c("character"))
  stopifnot( length(gene_id_genotypes) == length( snp_id_genotypes ) )
  
  stopifnot( class( sample_id ) %in% c("character"))
  
  ### keep genes that are in counts and in genotypes
  
  genes2keep <- gene_id_counts %in% gene_id_genotypes
  counts <- counts[genes2keep, , drop = FALSE]
  gene_id_counts <- gene_id_counts[genes2keep]
  feature_id_counts <- feature_id_counts[genes2keep]

  genes2keep <- gene_id_genotypes %in% gene_id_counts
  genotypes <- genotypes[genes2keep, , drop = FALSE]
  gene_id_genotypes <- gene_id_genotypes[genes2keep]
  snp_id_genotypes <- snp_id_genotypes[genes2keep]
  
  ### order genes in genotypes as in counts
  
  gene_id_genotypes <- factor(gene_id_genotypes, levels = unique(gene_id_counts))
  gene_id_counts <- factor(gene_id_counts, levels = unique(gene_id_counts))
  
  
  order_genotypes <- order(gene_id_genotypes)
  genotypes <- genotypes[order_genotypes, , drop = FALSE]
  gene_id_genotypes <- gene_id_genotypes[order_genotypes]
  snp_id_genotypes <- snp_id_genotypes[order_genotypes]
  
  colnames(counts) <- sample_id
  rownames(counts) <- feature_id_counts
  
  colnames(genotypes) <- sample_id
  rownames(genotypes) <- snp_id_genotypes
  
  counts <- new( "MatrixList", unlistData = counts, partitioning = IRanges::PartitioningByEnd(as.numeric(gene_id_counts), NG = nlevels(gene_id_counts), names = levels(gene_id_counts)) )
  
  genotypes <- new( "MatrixList", unlistData = genotypes, partitioning = IRanges::PartitioningByEnd(as.numeric(gene_id_genotypes), NG = nlevels(gene_id_genotypes), names = levels(gene_id_genotypes)) )
  
  samples <- IRanges::DataFrame(sample_id = sample_id)
  
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

  genotypes <- genotypes[s, ]
  snp_id_genotypes <- snp_id_genotypes[s]
  gene_id_genotypes <- names(gene_ranges)[q]


  data <- dmSQTLdata(counts = counts, gene_id_counts = gene_id_counts, feature_id_counts = feature_id_counts, genotypes = genotypes, gene_id_genotypes = gene_id_genotypes, snp_id_genotypes = snp_id_genotypes, sample_id = sample_id)

  return(data)

}


setMethod("show", "dmSQTLdata", function(object){
  
  cat("An object of class", class(object), "\n")
  
  cat("Slot \"counts\":\n")
  print(object@counts)
  
  cat("\nSlot \"genotypes\":\n")
  print(object@genotypes)
  
  cat("\nSlot \"samples\":\n")
  print(object@samples)
  
  })




setMethod("names", "dmSQTLdata", function(x) names(data@counts) )


setGeneric("dmSQTLfilter", function(x, ...) standardGeneric("dmSQTLfilter"))


setMethod("dmSQTLfilter", "dmSQTLdata", function(x, min_samps_gene_expr = 70, min_gene_expr = 1, min_samps_feature_prop = 5, min_feature_prop = 0.1, max_features = Inf, minor_allel_freq = 0.05, BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  # counts = x@counts; genotypes = x@genotypes; samples = x@samples
  
  data_filtered <- dmSQTL_filter(counts = x@counts, genotypes = x@genotypes, samples = x@samples, min_samps_gene_expr = min_samps_gene_expr, min_gene_expr = min_gene_expr, min_samps_feature_prop = min_samps_feature_prop, min_feature_prop = min_feature_prop, max_features = max_features, minor_allel_freq = minor_allel_freq, BPPARAM = BPPARAM)
  
  return(data_filtered)
  
  
  })



setGeneric("dmSQTLplotData", function(x, ...) standardGeneric("dmSQTLplotData"))

setMethod("dmSQTLplotData", "dmSQTLdata", function(x, out_dir = NULL){
  
 dmSQTL_plotData(counts = x@counts, genotypes = x@genotypes, out_dir = out_dir)
  
  })































