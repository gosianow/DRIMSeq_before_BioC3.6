#' @include class_MatrixList.R
NULL

##############################################################

#' Object that contains counts, genotypes and sample information.
#' 
#' Can be created with function \code{\link{dmSQTLdata}} and \code{\link{dmSQTLdataFromRanges}}.
#' 
#' @slot counts \code{\linkS4class{MatrixList}}  of counts.
#' @slot genotypes \code{\linkS4class{MatrixList}} with genotypes.
#' @slot samples data.frame with information about samples. Contains unique sample names (\code{sample_id}). 
setClass("dmSQTLdata", 
         representation(counts = "MatrixList", 
                        genotypes = "MatrixList", 
                        samples = "data.frame"))


##############################################################



##############################################################

setMethod("show", "dmSQTLdata", function(object){
  
  cat("An object of class", class(object), "\n")
  
  cat("with", length(object), "genes and", ncol(object@counts), "samples\n")
  
})



##############################################################

#' @rdname dmSQTLdata-class
#' @export
setMethod("names", "dmSQTLdata", function(x) names(x@counts) )

#' @rdname dmSQTLdata-class
#' @export
setMethod("length", "dmSQTLdata", function(x) length(x@counts) )


#' @rdname dmSQTLdata-class
#' @export
setMethod("[", "dmSQTLdata", function(x, i, j){
  
  if(missing(j)){
    
    counts <- x@counts[i, ]
    genotypes <- x@genotypes[i, ]
    samples <- x@samples
    
  }else{
    
    counts <- x@counts[i, j]
    genotypes <- x@genotypes[i, j]
    samples <- x@samples
    rownames(samples) <- samples$sample_id
    samples <- samples[j, ]
    samples$sample_id <- factor(samples$sample_id)
    rownames(samples) <- NULL
    
  }

  return(new("dmSQTLdata", counts = counts, genotypes = genotypes, samples = samples))
  
})

##############################################################

#'  Create \code{\linkS4class{dmSQTLdata}} object from tables of counts and matched genotypes
#'  
#'  @inheritParams dmDSdata
#'  @param genotypes A numeric matrix of matched genotypes. See Details.
#'  @param gene_id_genotypes Vector of gene IDs that are matched with genotypes.
#'  @param snp_id Vector of SNP IDs that correspond to genotypes.
#'  @return Returns a \code{\linkS4class{dmSQTLdata}} object containing counts, genotypes and sample information.
#'  @examples 
#'  ### counts
#'  head(dataSQTL_counts)
#'  
#'  counts <- dataSQTL_counts[, -(1:2)]
#'  gene_id <- group_split[, 1]
#'  feature_id <- group_split[, 2]
#'  
#'  ### gene_ranges
#'  dataSQTL_gene_ranges
#'  gene_ranges <- dataSQTL_gene_ranges
#'  names(gene_ranges) <- S4Vectors::mcols(gene_ranges)$name
#'  
#'  ### genotypes
#'  head(dataSQTL_genotypes)
#'  genotypes <- dataSQTL_genotypes[, -(1:4)]
#'  
#'  snp_id <- dataSQTL_genotypes$snp_id
#'  
#'  snp_ranges <- GenomicRanges::GRanges(S4Vectors::Rle(dataSQTL_genotypes$chr), IRanges::IRanges(dataSQTL_genotypes$start, dataSQTL_genotypes$end))
#'  names(snp_ranges) <- dataSQTL_genotypes$snp_id
#'  
#'  all(colnames(counts) == colnames(genotypes))
#'  
#'  sample_id <- colnames(counts)
#'  
#'  ### create dmSQTLdata object 
#'  d <- dmSQTLdataFromRanges(counts, gene_id, feature_id, gene_ranges, genotypes, snp_id, snp_ranges, sample_id, window = 5e3)
#'  
#'  @export
dmSQTLdata <- function(counts, gene_id, feature_id, genotypes, gene_id_genotypes, snp_id, sample_id){
  
  stopifnot( class( counts ) %in% c("matrix", "data.frame"))
  counts <- as.matrix(counts)
  stopifnot( mode( counts ) %in% c("numeric"))
  # counts <- ceiling(counts)
  
  stopifnot( class( gene_id ) %in% c("character", "factor"))
  stopifnot( class( feature_id ) %in% c("character", "factor"))
  stopifnot( length(gene_id) == length( feature_id ) )
  
  stopifnot( class( genotypes ) %in% c("matrix", "data.frame"))
  genotypes <- as.matrix(genotypes)
  stopifnot( mode( genotypes ) %in% c("numeric"))
  
  stopifnot( class( gene_id_genotypes ) %in% c("character", "factor"))
  stopifnot( class( snp_id ) %in% c("character", "factor"))
  stopifnot( length(gene_id_genotypes) == length( snp_id ) )
  
  stopifnot( class( sample_id ) %in% c("character", "factor"))
  
  ### keep genes that are in counts and in genotypes
  
  genes2keep <- gene_id %in% gene_id_genotypes
  counts <- counts[genes2keep, , drop = FALSE]
  gene_id <- gene_id[genes2keep]
  feature_id <- feature_id[genes2keep]
  
  genes2keep <- gene_id_genotypes %in% gene_id
  genotypes <- genotypes[genes2keep, , drop = FALSE]
  gene_id_genotypes <- gene_id_genotypes[genes2keep]
  snp_id <- snp_id[genes2keep]
  
  
  
  ### order genes in counts and in genotypes
  if(class(gene_id) == "character")
  gene_id <- factor(gene_id, levels = unique(gene_id))
  
  order_counts <- order(gene_id)
  counts <- counts[order_counts, , drop = FALSE]
  gene_id <- gene_id[order_counts]
  feature_id <- feature_id[order_counts]
  
  
  gene_id_genotypes <- factor(gene_id_genotypes, levels = levels(gene_id))

  order_genotypes <- order(gene_id_genotypes)
  genotypes <- genotypes[order_genotypes, , drop = FALSE]
  gene_id_genotypes <- gene_id_genotypes[order_genotypes]
  snp_id <- snp_id[order_genotypes]
  
  
  colnames(counts) <- sample_id
  rownames(counts) <- feature_id
  
  colnames(genotypes) <- sample_id
  rownames(genotypes) <- snp_id
  
  
  inds_counts <- 1:length(gene_id)
  names(inds_counts) <- feature_id
  partitioning_counts <- split(inds_counts, gene_id)
  
  inds_genotypes <- 1:length(gene_id_genotypes)
  names(inds_genotypes) <- snp_id
  partitioning_genotypes <- split(inds_genotypes, gene_id_genotypes)
  
  
  counts <- new( "MatrixList", unlistData = counts, partitioning = partitioning_counts)
  
  genotypes <- new( "MatrixList", unlistData = genotypes, partitioning = partitioning_genotypes)
  
  samples <- data.frame(sample_id = sample_id)
  
  data <- new("dmSQTLdata", counts = counts, genotypes = genotypes, samples = samples)
  
  return(data)
  
}

##############################################################

#' @rdname dmSQTLdata
#'  @param gene_ranges \code{\linkS4class{GRanges}} object with information 
#'    about gene location.
#'  @param genotypes A numeric matrix of unmatched genotypes. See Details.
#'  @param snp_ranges \code{\linkS4class{GRanges}} object with information about
#'    SNP location.
#'  @param window Numeric. Size of a down and up stream window that is used to 
#'    match SNPs to a gene. See details.
#'  @export
dmSQTLdataFromRanges <- function(counts, gene_id, feature_id, gene_ranges, genotypes, snp_id, snp_ranges, sample_id, window = 5e3){
  
  rownames(genotypes) <- snp_id
  gene_ranges <- GenomicRanges::resize(gene_ranges, GenomicRanges::width(gene_ranges) + 2 * window, fix = "center")
  
  ## Match genes and SNPs
  variantMatch <- GenomicRanges::findOverlaps(gene_ranges, snp_ranges, select = "all")
  
  q <- GenomicRanges::queryHits(variantMatch)
  s <- GenomicRanges::subjectHits(variantMatch)
  
  genotypes <- genotypes[s, ]
  snp_id <- snp_id[s]
  gene_id_genotypes <- names(gene_ranges)[q]
  
  
  data <- dmSQTLdata(counts = counts, gene_id = gene_id, feature_id = feature_id, genotypes = genotypes, gene_id_genotypes = gene_id_genotypes, snp_id = snp_id, sample_id = sample_id)
  
  return(data)
  
}



##############################################################

#' @rdname dmFilter
#' @inheritParams dmFilter
#' @param minor_allel_freq Number of samples that the minor allel must be present in.
#' @param BPPARAM Parallelization method used by
#'   \code{\link[BiocParallel]{bplapply}}.
#' @export
setMethod("dmFilter", "dmSQTLdata", function(x, min_samps_gene_expr = 70, min_gene_expr = 1, min_samps_feature_prop = 5, min_feature_prop = 0.1, max_features = Inf, minor_allel_freq = 5, BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  data_filtered <- dmSQTL_filter(counts = x@counts, genotypes = x@genotypes, samples = x@samples, min_samps_gene_expr = min_samps_gene_expr, min_gene_expr = min_gene_expr, min_samps_feature_prop = min_samps_feature_prop, min_feature_prop = min_feature_prop, max_features = max_features, minor_allel_freq = minor_allel_freq, BPPARAM = BPPARAM)
  
  return(data_filtered)
  
  
})



##############################################################

#' @rdname plotData
#' @export
setMethod("plotData", "dmSQTLdata", function(x, out_dir = NULL){
  
  dmSQTL_plotData(counts = x@counts, genotypes = x@genotypes, out_dir = out_dir)
  
})



##############################################################

#' @rdname dmSQTLdata-class
#' @export
setMethod("plot", "dmSQTLdata", function(x, out_dir = NULL){
  
  plotData(x, out_dir = out_dir)
  
})




























