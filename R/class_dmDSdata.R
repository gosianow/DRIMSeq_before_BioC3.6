#' @include class_MatrixList.R
NULL


##############################################################
#' dmDSdata object
#' 
#' dmDSdata contains expression, in counts, of genomic features such as exons or transcripts and sample information needed for the differential splicing (DS) analysis. It can be created with function \code{\link{dmDSdata}}.
#' 
#' @details 
#' 
#' \itemize{
#'  \item \code{counts(object)}: Get a data.frame with counts.
#'  \item \code{samples(x)}: Get a data.frame with the sample information.
#'   \item \code{names(x)}: Get the gene names.
#'   \item \code{length(x)}: Get the number of genes.
#'   \item \code{x[i, j]}: Get a subset of dmDSdata object that consists of counts for genes i and samples j. 
#' }
#' 
#' @param object,x dmDSdata object.
#' @param i,j Parameters used for subsetting. See Details.
#' @param ... Other parameters that can be defined by methods using this generic.
#' 
#' @slot counts \code{\linkS4class{MatrixList}} of expression, in counts, of genomic features. Rows correspond to genomic features, such as exons or transcripts. Columns correspond to samples. MatrixList is partitioned in a way that each of the matrices in a list contains counts for a single gene.
#' @slot samples data.frame with information about samples. It must contain variables: \code{sample_id} of unique sample names and \code{group} which groups samples into conditions.
#' @author Malgorzata Nowicka
#' @seealso \code{\link{plotData}}, \code{\linkS4class{dmDSdispersion}}, \code{\linkS4class{dmDSfit}}, \code{\linkS4class{dmDStest}}
setClass("dmDSdata", 
         representation(counts = "MatrixList", 
                        samples = "data.frame"))



setValidity("dmDSdata", function(object){
  # has to return TRUE when valid object!
  
  if(ncol(object@counts) == nrow(object@samples))
    out <- TRUE
  else 
    return(paste0("Unequal number of samples in 'counts' and 'samples' ", ncol(object@counts), " and ", nrow(object@samples)))
  
  if(all(c("sample_id", "group") %in% colnames(object@samples)))
    out <- TRUE
  else
    return(paste0("'samples' must contain 'sample_id' amd 'group' variables"))
  
  if(length(unique(object@samples$sample_id)) == nrow(object@samples))
    out <- TRUE
  else
    return("There must be a unique 'sample_id' for each sample")
  
  return(out)
  
})


##############################################################


#' @rdname dmDSdata-class
#' @export
setMethod("counts", "dmDSdata", function(object){
  
  data.frame(gene_id = rep(names(object@counts), width(object@counts)), feature_id = rownames(object@counts@unlistData), object@counts@unlistData, stringsAsFactors = FALSE, row.names = NULL)
  
})


#' @rdname dmDSdata-class
#' @export
setGeneric("samples", function(x, ...) standardGeneric("samples"))


#' @rdname dmDSdata-class
#' @export
setMethod("samples", "dmDSdata", function(x) x@samples )



##############################################################

setMethod("show", "dmDSdata", function(object){
  
  cat("An object of class", class(object), "\n")
  
  cat("with", length(object), "genes and", ncol(object@counts), "samples\n")
  
  cat("* data accessors: counts(), samples()\n")
  
})

##############################################################

#' @rdname dmDSdata-class
#' @export
setMethod("names", "dmDSdata", function(x) names(x@counts) )

#' @rdname dmDSdata-class
#' @export
setMethod("length", "dmDSdata", function(x) length(x@counts) )


#' @aliases [,dmDSdata-method
#' @rdname dmDSdata-class
#' @export
setMethod("[", "dmDSdata", function(x, i, j){
  
  if(missing(j)){
    
    counts <- x@counts[i, ]
    samples <- x@samples
    
  }else{
    
    counts <- x@counts[i, j]
    samples <- x@samples
    rownames(samples) <- samples$sample_id
    samples <- samples[j, ]
    samples$sample_id <- factor(samples$sample_id)
    samples$group <- factor(samples$group)
    rownames(samples) <- NULL
    
  }
  
  return(new("dmDSdata", counts = counts, samples = samples))
  
})



##############################################################

#'  Create dmDSdata object
#'  
#'  Constructor function for a \code{\linkS4class{dmDSdata}} object. 
#'  
#'  @param counts Numeric matrix or data.frame of counts. Rows represent features, for example, exons,
#'    exonic bins or transcripts. Columns represent samples.
#'  @param gene_id Vector of gene IDs corresponding to \code{counts}.
#'  @param feature_id Vector of feature IDs corresponding to \code{counts}.
#'  @param sample_id Vector of unique sample IDs corresponding to the columns in \code{counts}.
#'  @param group Vector that defines the goupping of samples.
#'  @return Returns a \linkS4class{dmDSdata} object.
#'  @examples
#'  ### Differential splicing analysis
#'  
#'  head(dataDS_counts)
#'  dataDS_metadata
#'  
#'  counts <- dataDS_counts[,-1]
#'  group_id <- dataDS_counts[,1]
#'  group_split <- limma::strsplit2(group_id, ":")
#'  gene_id <- group_split[, 1]
#'  feature_id <- group_split[, 2]
#'  sample_id = dataDS_metadata$sample_id
#'  group = dataDS_metadata$group
#'  
#'  d <- dmDSdata(counts = counts, gene_id = gene_id, feature_id = feature_id, 
#'  sample_id = sample_id, group = group)
#'  
#'  length(d)
#'  
#'  plotData(d)
#'  
#'  @seealso \code{\link{plotData}}, \code{\link{dmFilter}}, \code{\link{dmDispersion}}, \code{\link{dmFit}}, \code{\link{dmTest}}
#'  @author Malgorzata Nowicka
#'  @export
dmDSdata <- function(counts, gene_id, feature_id, sample_id, group){
  
  stopifnot(class(counts) %in% c("matrix", "data.frame"))
  counts <- as.matrix(counts)
  stopifnot(mode(counts) %in% "numeric")
  # counts <- ceiling(counts)
  
  stopifnot( class( gene_id ) %in% c("character", "factor"))
  stopifnot( class( feature_id ) %in% c("character", "factor"))
  stopifnot( class( sample_id ) %in% c("character", "factor"))
  stopifnot( class( group ) %in% c("character", "factor"))
  stopifnot( length(gene_id) == length( feature_id ) )
  stopifnot( length(sample_id) == length( group ) )
  stopifnot(nrow(counts) == length(feature_id))
  stopifnot(ncol(counts) == length(sample_id))
  
  if(class(gene_id) == "character")
    gene_id <- factor(gene_id, levels = unique(gene_id))
  else 
    gene_id <- factor(gene_id)
  
  if(class(group) == "character")
    group <- factor(group, levels = unique(group))
  else 
    group <- factor(group)
  
  if(class(sample_id) == "character")
    sample_id <- factor(sample_id, levels = unique(sample_id))
  else 
    sample_id <- factor(sample_id)
  
  
  tbl <- table(group)
  if(all(tbl < 2))
    stop("There must be at least one group with two replicates!")
  
  ### ordering
  or <- order(gene_id)
  oc <- order(group)
  
  counts <- counts[or, oc, drop = FALSE]
  gene_id <- gene_id[or]
  feature_id <- feature_id[or]
  group <- group[oc]
  sample_id <- sample_id[oc]
  
  colnames(counts) <- sample_id
  rownames(counts) <- feature_id
  
  inds <- 1:length(gene_id)
  names(inds) <- feature_id
  
  partitioning <- split(inds, gene_id)
  
  samples <- data.frame(sample_id = sample_id, group = group)
  
  data <- new("dmDSdata", counts = new("MatrixList", unlistData = counts, partitioning = partitioning), samples = samples)
  
  return(data)
  
}


##############################################################

#' Filtering
#' 
#' Filtering of genes and features with too low expression. Additionally, for the dmSQTLdata object, filtering of genotypes with too low frequency.
#' 
#' @param x \code{\linkS4class{dmDSdata}} or \code{\linkS4class{dmSQTLdata}} object.
#' @param ... Other parameters that can be defined by methods using this generic.
#' @export
setGeneric("dmFilter", function(x, ...) standardGeneric("dmFilter"))



##############################################################

#' @details 
#' Filtering parameters should be adjusted according to the sample size of the experiment data and the number of replicates per condition. 
#' 
#' In differential splicing analysis, we suggest using \code{min_samps_gene_expr} and \code{min_samps_feature_prop} equal to the mininal number of replicates in any of the conditions. For example, in an assay with 3 versus 5 replicates, we would set this parameters to 3.
#' 
#' @param min_samps_gene_expr Minimal number of samples where genes should 
#'   be expressed. See Details.
#' @param min_gene_expr Minimal gene expression in CPM \code{\link[edgeR]{cpm}}.
#' @param min_samps_feature_prop Minimal number of samples where features
#'   should be expressed. See details.
#' @param min_feature_prop Minimal proportion for feature expression. This value
#'   should be between 0 and 1.
#' @param max_features Maximum number of features, which pass the filtering criteria, that should be kept for each gene. If equal to \code{Inf}, all features that pass the filtering criteria are kept. 
#' @return Returns filtered \code{\linkS4class{dmDSdata}} or \code{\linkS4class{dmSQTLdata}} object.
#' @examples 
#' ### Differential splicing analysis
#'
#' dd <- dataDS_dmDSdata
#' plotData(dd)
#' 
#' d <- dmFilter(dd)
#' plotData(d)
#' 
#' d <- dmFilter(dd, max_features = 10)
#' plotData(d)
#' 
#' @seealso \code{\link{plotData}}, \code{\link{dmDispersion}}, \code{\link{dmFit}}, \code{\link{dmTest}}
#' @author Malgorzata Nowicka
#' @rdname dmFilter
#' @export
setMethod("dmFilter", "dmDSdata", function(x, min_samps_gene_expr = 3, min_gene_expr = 1, min_samps_feature_prop = 3, min_feature_prop = 0.01, max_features = Inf){
  
  data_filtered <- dmDS_filter(counts = x@counts, samples = x@samples, min_samps_gene_expr = min_samps_gene_expr, min_gene_expr = min_gene_expr, min_samps_feature_prop = min_samps_feature_prop, min_feature_prop = min_feature_prop, max_features = max_features)
  
  return(data_filtered)
  
  
})

##############################################################

#' Plot data summary
#' 
#' Plot a histogram of the number of features per gene. Additionally, for \code{\linkS4class{dmSQTLdata}} object, plot a histogram of the number of SNPs per gene and a histogram of the number of unique SNPs (blocks) per gene.
#' 
#' @param x \code{\linkS4class{dmDSdata}} or \code{\linkS4class{dmSQTLdata}} object.
#' @param ... Other parameters that can be defined by methods using this generic.
#' @export
setGeneric("plotData", function(x, ...) standardGeneric("plotData"))


##############################################################

#' @param out_dir Directory where the plot should be saved. If \code{NULL} the plot is printed. 
#' @examples 
#' ### Differential splicing analysis
#'
#' d <- dataDS_dmDSdata
#' plotData(d)
#'
#' @author Malgorzata Nowicka
#' @seealso \code{\link{plotDispersion}}, \code{\link{plotFit}}, \code{\link{plotTest}}
#' @rdname plotData
#' @export
setMethod("plotData", "dmDSdata", function(x, out_dir = NULL){
  
  dmDS_plotData(counts = x@counts, out_dir = out_dir)
  
})



