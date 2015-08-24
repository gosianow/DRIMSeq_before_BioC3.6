#' @include class_MatrixList.R
NULL


##############################################################
#' Object that contains counts and sample information.
#' 
#' Can be created with function \code{\link{dmDSdata}}.
#' 
#' @slot counts \code{\linkS4class{MatrixList}} of counts.
#' @slot samples data.frame with information about samples. Contains unique sample names (\code{sample_id}) and information about grouping into conditions (\code{group}).
setClass("dmDSdata", 
         representation(counts = "MatrixList", samples = "data.frame"))


##############################################################
#'  Create dmDSdata object from a table of counts
#'  
#'  @param counts A numeric matrix or data.frame of counts. Rows represent features (exons,
#'    bins or transcripts), columns represent samples.
#'  @param gene_id_counts Vector of gene IDs of lenght correspoding to the number of rows in \code{counts}.
#'  @param feature_id_counts Vector of feature IDs of lenght correspoding to the number of rows in \code{counts}.
#'  @param sample_id A vector of unique sample IDs of length corresponding to
#'    the number of columns in \code{counts}.
#'  @param group A vector that defines the goupping of samples.
#'  @return Returns \code{dmDSdata} object containing counts and sample 
#'    information \code{\linkS4class{dmDSdata}}.
#'  @examples
#'  head(dataDS_counts)
#'  dataDS_metadata
#'  
#'  counts <- as.matrix(dataDS_counts[,-1])
#'  group_id <- dataDS_counts[,1]
#'  group_split <- limma::strsplit2(group_id, ":")
#'  gene_id_counts <- group_split[, 1]
#'  feature_id_counts <- group_split[, 2]
#'  sample_id = dataDS_metadata$sample_id
#'  group = dataDS_metadata$group
#'  
#'  data <- dmDSdata(counts = counts, gene_id_counts = gene_id_counts, feature_id_counts = feature_id_counts, sample_id = sample_id, group = group)
#'  dmDSplotData(data)
#'  @export
dmDSdata <- function(counts, gene_id_counts, feature_id_counts, sample_id, group){
  
  stopifnot(class(counts) %in% c("matrix", "data.frame"))
  counts <- as.matrix(counts)
  stopifnot(mode(counts) %in% "numeric")
  counts <- ceiling(counts)
  
  stopifnot( class( gene_id_counts ) %in% c("character", "factor"))
  stopifnot( class( feature_id_counts ) %in% c("character", "factor"))
  stopifnot( class( sample_id ) %in% c("character", "factor"))
  stopifnot( class( group ) %in% c("character", "factor"))
  stopifnot( length(gene_id_counts) == length( feature_id_counts ) )
  stopifnot( length(sample_id) == length( group ) )
  stopifnot(nrow(counts) == length(feature_id_counts))
  stopifnot(ncol(counts) == length(sample_id))
  
  if(class(gene_id_counts) == "character")
  gene_id_counts <- factor(gene_id_counts, levels = unique(gene_id_counts))
  
  if(class(group) == "character")
  group <- factor(group, levels = unique(group))
  
  
  ### keep samples/groups with enough replicates
  tbl <- table(group)
  
  if(sum(tbl > 1) < 2)
  stop("There must be at least two groups and every group with at least two replicates!")
  
  levels <- names(tbl[tbl > 1])
  
  if(length(levels) < length(tbl))
  message("Groups that are kept: ", paste(levels, collapse = ", "))
  
  keep_samps <- group %in% levels
  
  group <- factor(group[keep_samps], levels = levels)
  
  counts <- counts[, keep_samps]
  sample_id <- sample_id[keep_samps]
  
  
  ### ordering
  or <- order(gene_id_counts)
  oc <- order(group)
  
  counts <- counts[or, oc, drop = FALSE]
  gene_id_counts <- gene_id_counts[or]
  feature_id_counts <- feature_id_counts[or]
  group <- group[oc]
  sample_id <- sample_id[oc]
  
  colnames(counts) <- sample_id
  rownames(counts) <- feature_id_counts
  
  partitioning <- split(1:length(gene_id_counts), gene_id_counts)

  samples <- data.frame(sample_id = sample_id, group = group)
  
  data <- new("dmDSdata", counts = new("MatrixList", unlistData = counts, partitioning = partitioning), samples = samples)
  
  return(data)
  
}

# setValidity("dmDSdata", function(object){

#   # has to return TRUE when valid object!

#   })

##############################################################
setMethod("show", "dmDSdata", function(object){
  
  cat("An object of class", class(object), "\n")
  
  cat("Slot \"counts\":\n")
  print(object@counts)
  
  cat("\nSlot \"samples\":\n")
  show_matrix(object@samples)
  
})

##############################################################
setMethod("names", "dmDSdata", function(x) names(x@counts) )


##############################################################
#' Filtering.
#' 
#' Filtering of genes with low expression and features with low proportions.
#' 
#' Parameters should be adjusted based on the sample size used for the analysis...
#' 
#' @param x \code{\link{dmDSdata}} object of counts.
#' @param ... Filtering parameters.
#' @export
setGeneric("dmDSfilter", function(x, ...) standardGeneric("dmDSfilter"))

##############################################################
#' @rdname dmDSfilter
#' @param min_samps_gene_expr Minimal number of samples where the genes should 
#'   be expressed.
#' @param min_gene_expr Minimal gene expression in CPM \code{\link[edgeR]{cpm}}.
#' @param min_samps_feature_prop Minimal number of samples where the features
#'   should be expressed.
#' @param min_feature_prop Minimal proportion for feature expression. This value
#'   should be between 0 and 1.
#' @param max_features Maximum number of features that should be kept per gene.
#' @return Returns filtered \code{\linkS4class{dmDSdata}} object.
#' @examples 
#' data <- dataDS_dmDSdata
#' dmDSplotData(data)
#' data <- dmDSfilter(data)
#' dmDSplotData(data)
#' @export
setMethod("dmDSfilter", "dmDSdata", function(x, min_samps_gene_expr = 3, min_gene_expr = 1, min_samps_feature_prop = 3, min_feature_prop = 0.01, max_features = Inf){
  
  
  data_filtered <- dmDS_filter(counts = x@counts, samples = x@samples, min_samps_gene_expr = min_samps_gene_expr, min_gene_expr = min_gene_expr, min_samps_feature_prop = min_samps_feature_prop, min_feature_prop = min_feature_prop, max_features = max_features)
  
  return(data_filtered)
  
  
})

##############################################################
#' Plot the data information.
#' 
#' Plots a histogram of number of features per gene.
#' 
#' @param x \code{\link{dmDSdata}} object of counts.
#' @param ... Plotting parameters.
#' @export
setGeneric("dmDSplotData", function(x, ...) standardGeneric("dmDSplotData"))


##############################################################
#' @rdname dmDSplotData
#' @param out_dir Directory where the plot should be saved. If \code{NULL} the plot is printed.
#' @param info \code{DataFrame} with \code{gene_id} and \code{feature_id} that are differentially spliced (DS). 
#' @examples 
#' dmDSplotData(dataDS_dmDSdata)
#' @export
setMethod("dmDSplotData", "dmDSdata", function(x, out_dir = NULL, info = NULL){
  
  dmDS_plotData(counts = x@counts, out_dir = out_dir, info = info)
  
})


















































