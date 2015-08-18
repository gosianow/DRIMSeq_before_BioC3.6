#' @include class_MatrixList.R
NULL

##############################################################
#' Object that contains counts, genotypes and sample information.
#' 
#' Can be created with function \code{\link{dmSQTLdata}} and \code{\link{dmSQTLdataFromRanges}}.
#' 
#' @slot counts \code{\linkS4class{MatrixList}}  of counts.
#' @slot genotypes \code{\linkS4class{MatrixList}} with genotypes.
#' @slot samples DataFrame with information about samples. Contains unique sample names (\code{sample_id}). 
#' @importClassesFrom S4Vectors DataFrame
setClass("dmSQTLdata", 
         representation(counts = "MatrixList", 
                        genotypes = "MatrixList", 
                        samples = "DataFrame"))



##############################################################
#'  Create \code{\linkS4class{dmSQTLdata}} object from tables of counts and matched genotypes
#'  
#'  @inheritParams dmDSdata
#'  @param genotypes A numeric matrix of matched genotypes. See Details.
#'  @param gene_id_genotypes Vector of gene IDs that are matched with genotypes.
#'  @param snp_id_genotypes Vector of SNP IDs that correspond to genotypes.
#'  @return Returns a \code{\linkS4class{dmSQTLdata}} object containing counts, genotypes and sample 
#'    information.
#'  @export
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
  
  samples <- S4Vectors::DataFrame(sample_id = sample_id)
  
  data <- new("dmSQTLdata", counts = counts, genotypes = genotypes, samples = samples)
  
  return(data)
  
}

##############################################################

#'  Create \code{\linkS4class{dmSQTLdata}} object from tables of counts, 
#'  genotypes and gene ranges
#'  
#'  @inheritParams dmSQTLdata
#'  @param gene_ranges \code{\linkS4class{GRanges}} object with information 
#'    about gene location.
#'  @param genotypes A numeric matrix of unmatched genotypes. See Details.
#'  @param snp_ranges \code{\linkS4class{GRanges}} object with information about
#'    SNP location.
#'  @param window Numeric. Size of a down and up stream window that is used to 
#'    match SNPs to a gene. See details.
#'  @return Returns a \code{\linkS4class{dmSQTLdata}} object containing counts, 
#'    genotypes and sample information.
#'    @examples 
#'    ### counts
#'  head(dataSQTL_counts)
#'  counts <- as.matrix(dataSQTL_counts[, -1])
#'  
#'  group_id <- dataSQTL_counts[,1]
#'  group_split <- limma::strsplit2(group_id, ":")
#'  gene_id_counts <- group_split[, 1]
#'  feature_id_counts <- group_split[, 2]
#'  
#'  ### gene_ranges
#'  dataSQTL_gene_ranges
#'  gene_ranges <- dataSQTL_gene_ranges
#'  names(gene_ranges) <- S4Vectors::mcols(gene_ranges)$name
#'  
#'  
#'  ### genotypes
#'  head(dataSQTL_genotypes)
#'  genotypes <- as.matrix(dataSQTL_genotypes[, -(1:4)])
#'  
#'  snp_id_genotypes <- dataSQTL_genotypes$snp_id
#'  
#'  snp_ranges <- GenomicRanges::GRanges(S4Vectors::Rle(dataSQTL_genotypes$chr), 
#'  IRanges::IRanges(dataSQTL_genotypes$start, dataSQTL_genotypes$end))
#'  names(snp_ranges) <- dataSQTL_genotypes$snp_id
#'  
#'  all(colnames(counts) == colnames(genotypes))
#'  
#'  sample_id <- colnames(counts)
#'  
#'  ### create dmSQTLdata object 
#'  data <- dmSQTLdataFromRanges(counts, gene_id_counts, feature_id_counts, 
#'  gene_ranges, genotypes, snp_id_genotypes, snp_ranges, sample_id, 
#'  window = 5e3)
#'  
#'  dmSQTLplotData(data)
#'  
#'  @export
dmSQTLdataFromRanges <- function(counts, gene_id_counts, feature_id_counts, gene_ranges, genotypes, snp_id_genotypes, snp_ranges, sample_id, window = 5e3){
  
  rownames(genotypes) <- snp_id_genotypes
  gene_ranges <- GenomicRanges::resize(gene_ranges, GenomicRanges::width(gene_ranges) + 2 * window, fix = "center")
  
  ## Match genes and SNPs
  variantMatch <- GenomicRanges::findOverlaps(gene_ranges, snp_ranges, select = "all")
  
  q <- GenomicRanges::queryHits(variantMatch)
  s <- GenomicRanges::subjectHits(variantMatch)
  
  genotypes <- genotypes[s, ]
  snp_id_genotypes <- snp_id_genotypes[s]
  gene_id_genotypes <- names(gene_ranges)[q]
  
  
  data <- dmSQTLdata(counts = counts, gene_id_counts = gene_id_counts, feature_id_counts = feature_id_counts, genotypes = genotypes, gene_id_genotypes = gene_id_genotypes, snp_id_genotypes = snp_id_genotypes, sample_id = sample_id)
  
  return(data)
  
}

##############################################################
setMethod("show", "dmSQTLdata", function(object){
  
  cat("An object of class", class(object), "\n")
  
  cat("Slot \"counts\":\n")
  print(object@counts)
  
  cat("\nSlot \"genotypes\":\n")
  print(object@genotypes)
  
  cat("\nSlot \"samples\":\n")
  print(object@samples)
  
})



##############################################################
setMethod("names", "dmSQTLdata", function(x) names(data@counts) )




##############################################################
#' Filtering.
#' 
#' Two step filtering. First, for counts - filtering of genes with low expression and features with low proportions. Second, for genotypes - filtering genotypes with too low monor allel frequency.
#' 
#' Parameters should be adjusted based on the sample size used for the analysis...
#' 
#' @param x \code{\link{dmSQTLdata}} object with counts and genotypes.
#' @param ... Filtering parameters.
#' @export
setGeneric("dmSQTLfilter", function(x, ...) standardGeneric("dmSQTLfilter"))




##############################################################

#' @rdname dmSQTLfilter
#' @inheritParams dmDSfilter
#' @param minor_allel_freq Value between 0 and 1 which corresponds to minimal
#'   allel frequency.
#' @param BPPARAM Parallelization method used by
#'   \code{\link[BiocParallel]{bplapply}}.
#' @return Returns filtered \code{\linkS4class{dmSQTLdata}} object.
#' @examples 
#' data <- dataSQTL_dmSQTLdata
#' dmSQTLplotData(data)
#' 
#' data <- dmSQTLfilter(data)
#' dmSQTLplotData(data)
#' @export
setMethod("dmSQTLfilter", "dmSQTLdata", function(x, min_samps_gene_expr = 70, min_gene_expr = 1, min_samps_feature_prop = 5, min_feature_prop = 0.1, max_features = Inf, minor_allel_freq = 0.05, BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  # counts = x@counts; genotypes = x@genotypes; samples = x@samples
  
  data_filtered <- dmSQTL_filter(counts = x@counts, genotypes = x@genotypes, samples = x@samples, min_samps_gene_expr = min_samps_gene_expr, min_gene_expr = min_gene_expr, min_samps_feature_prop = min_samps_feature_prop, min_feature_prop = min_feature_prop, max_features = max_features, minor_allel_freq = minor_allel_freq, BPPARAM = BPPARAM)
  
  return(data_filtered)
  
  
})


##############################################################
#' Plot the data information.
#' 
#' Plots histograms of number of features per gene and number of SNPs per gene.
#' 
#' @param x \code{\link{dmSQTLdata}} object of counts.
#' @param ... Plotting parameters.
#' @export
setGeneric("dmSQTLplotData", function(x, ...) standardGeneric("dmSQTLplotData"))



##############################################################
#' @rdname dmSQTLplotData
#' @param out_dir Directory where the plot should be saved. If \code{NULL} the plot is printed.
#' @export
setMethod("dmSQTLplotData", "dmSQTLdata", function(x, out_dir = NULL){
  
  dmSQTL_plotData(counts = x@counts, genotypes = x@genotypes, out_dir = out_dir)
  
})































