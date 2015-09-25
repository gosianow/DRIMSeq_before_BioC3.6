#' @include class_MatrixList.R
NULL


##############################################################
#' dmDSdata object
#' 
#' dmDSdata contains counts of genomic features, such as exons or transcrips, and sample information needed for differential splicing (DS) analysis. It can be created with function \code{\link{dmDSdata}}.
#' 
#' @slot counts A \linkS4class{MatrixList} of counts where rows correspond to genomic features, such as exons or transcripts, that are quantified. Columns correspond to samples. MatrixList is partitioned in a way that each submatrix contains counts for a single gene.
#' @slot samples A data.frame with information about samples. It must contain columns with unique sample names (\code{sample_id}) and information about grouping into conditions (\code{group}).

setClass("dmDSdata", 
         representation(counts = "MatrixList", 
                        samples = "data.frame"))



setValidity("dmDSdata", function(object){
  # has to return TRUE when valid object!
  
  if(ncol(object@counts) == nrow(object@samples))
    out <- TRUE
  else 
    return(paste0("Unequal number of samples in 'counts' and 'samples' ", ncol(object@counts), " and ", nrow(object@samples)))
  
  if(c("sample_id", "group") %in% colnames(object@samples))
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
#' @param object dmDSdata object.
#' @export
setMethod("counts", "dmDSdata", function(object, ...){
  
  data.frame(gene_id = rep(names(object@counts), width(object@counts)), feature_id = rownames(object@counts@unlistData), object@counts@unlistData, stringsAsFactors = FALSE, row.names = NULL)
  
})


#' @rdname dmDSdata-class
#' @param x Object to extract the \code{samples} slot from.
#' @param ... Other parameters that can be used for the extraction of desired attributes.
#' @export
setGeneric("samples", function(x, ...) standardGeneric("samples"))

#' @rdname dmDSdata-class
#' @param x dmDSdata object.
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


#' @rdname dmDSdata-class
#' @param i,j Indexes defining which genes (i) and which samples (j) should be subsetted. 
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

#'  Create dmDSdata object from a table of counts
#'  
#'  @param counts A numeric matrix or data.frame of counts. Rows represent features (exons,
#'    bins or transcripts), columns represent samples.
#'  @param gene_id Vector of gene IDs of lenght correspoding to the number of rows in \code{counts}.
#'  @param feature_id Vector of feature IDs of lenght correspoding to the number of rows in \code{counts}.
#'  @param sample_id A vector of unique sample IDs of length corresponding to
#'    the number of columns in \code{counts}.
#'  @param group A vector that defines the goupping of samples.
#'  @return Returns \code{dmDSdata} object containing counts and sample 
#'    information \code{\linkS4class{dmDSdata}}.
#'  @examples
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
#'  d <- dmDSdata(counts = counts, gene_id = gene_id, feature_id = feature_id, sample_id = sample_id, group = group)
#'  
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

#' Filtering.
#' 
#' Filtering of genes with low expression and features with low proportions.
#' 
#' Parameters should be adjusted based on the sample size used for the analysis...
#' 
#' @param x \code{\linkS4class{dmDSdata}} or \code{\linkS4class{dmSQTLdata}} object.
#' @param ... Filtering parameters.
#' @export
setGeneric("dmFilter", function(x, ...) standardGeneric("dmFilter"))



##############################################################

#' @rdname dmFilter
#' @param min_samps_gene_expr Minimal number of samples where the genes should 
#'   be expressed.
#' @param min_gene_expr Minimal gene expression in CPM \code{\link[edgeR]{cpm}}.
#' @param min_samps_feature_prop Minimal number of samples where the features
#'   should be expressed.
#' @param min_feature_prop Minimal proportion for feature expression. This value
#'   should be between 0 and 1.
#' @param max_features Maximum number of features that should be kept per gene.
#' @return Returns filtered \code{\linkS4class{dmDSdata}} or \code{\linkS4class{dmSQTLdata}} object.
#' @examples 
#' ### dmDSdata
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
#' ### dmSQTLdata
#'
#' dd <- dataSQTL_dmSQTLdata
#' plotData(dd)
#' 
#' d <- dmFilter(dd)
#' plotData(d)
#' @export
setMethod("dmFilter", "dmDSdata", function(x, min_samps_gene_expr = 3, min_gene_expr = 1, min_samps_feature_prop = 3, min_feature_prop = 0.01, max_features = Inf){
  
  data_filtered <- dmDS_filter(counts = x@counts, samples = x@samples, min_samps_gene_expr = min_samps_gene_expr, min_gene_expr = min_gene_expr, min_samps_feature_prop = min_samps_feature_prop, min_feature_prop = min_feature_prop, max_features = max_features)
  
  return(data_filtered)
  
  
})

##############################################################

#' Plot the data information.
#' 
#' Plots a histogram of number of features per gene.
#' 
#' @param x \code{\linkS4class{dmDSdata}} or \code{\linkS4class{dmSQTLdata}} object with data.
#' @param ... Plotting parameters.
#' @export
setGeneric("plotData", function(x, ...) standardGeneric("plotData"))


##############################################################

#' @rdname plotData
#' @param out_dir Directory where the plot should be saved. If \code{NULL} the plot is printed. 
#' @examples 
#' ### dmDSdata
#'
#' d <- dataDS_dmDSdata
#' plotData(d)
#' plot(d)
#'
#' ### dmSQTLdata
#'
#' d <- dataSQTL_dmSQTLdata
#' plotData(d)
#' plot(d)
#' @export
setMethod("plotData", "dmDSdata", function(x, out_dir = NULL){
  
  dmDS_plotData(counts = x@counts, out_dir = out_dir)
  
})


##############################################################

#' @rdname dmDSdata-class
#' @export
setMethod("plot", "dmDSdata", function(x, out_dir = NULL){
  
  plotData(x, out_dir = out_dir)
  
})















































