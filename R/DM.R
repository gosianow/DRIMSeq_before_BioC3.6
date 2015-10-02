#' @importFrom methods setClass setGeneric setMethod initialize
#' @importFrom BiocGenerics width counts
#' @import ggplot2
#' @importFrom reshape2 melt
NULL





################################################################################
### Data documentation
################################################################################



#' Data for differential splicing analysis
#'
#' A subset (100 randomly selected genes) of exonic bin counts computed with DEXSeq package in a simulation study based on Drosophila genome. In this study data was generated to mimic an RNA-seq assay for 6 samples separated into 2 conditions (3 versus 3 samples). Differential transcript usage was induced by swapping, between the 2 conditions, the expression of two most abundant transcripts for 1000 genes. Our subset consists of 20 genes with differential splicing (status 1) and 80 without (status 0). For more details about this simulations, see the reference in Source paragraph.
#' 
#' @format 
#' \code{dataDS_counts} contains exonic bin counts. A data.frame with 7 variables:
#' \itemize{
#'   \item \code{group_id}: Exonic bin IDs.
#'   \item \code{sample_1, sample_2, ...}: Quantification of the exonic bins in 6 samples.
#'   }
#'   
#' \code{dataDS_metadata} describes the samples. A data.frame with 2 variables:
#' \itemize{
#'   \item \code{sample_id}: Sample IDs.
#'   \item \code{group}: Grouping into two conditions.
#' }
#' 
#' \code{dataDS_status} contains the information about the differential splicing status of genes. A data.frame with 2 variables:
#' \itemize{
#'   \item \code{gene_id}: Gene IDs.
#'   \item \code{status}: 1 for genes with differential splicing and 0 for genes without DS.
#' }
#' 
#' \code{dataDS_dmDSdata}, \code{dataDS_dmDSdispersion} and \code{dataDS_dmDStest} are the \code{DM} package objects created through the differential splicing analysis pipeline. 
#' 
#' @source 
#' Soneson C, Matthes KL, Nowicka M, Law CW, Robinson MD. Differential transcript usage from RNA-seq data : isoform pre-filtering improves performance of count-based methods. 2015.
#' 
#' \url{http://imlspenticton.uzh.ch/robinson_lab/splicing_comparison/}
"dataDS_counts"



#' @rdname dataDS_counts
"dataDS_metadata"

#' @rdname dataDS_counts
"dataDS_status"

#' @rdname dataDS_counts
"dataDS_dmDSdata"

#' @rdname dataDS_counts
"dataDS_dmDSdispersion"

#' @rdname dataDS_counts
"dataDS_dmDStest"







#' Data for sQTL analysis
#'
#' A subset of data from GEUVADIS project where 465 RNA-seq samples from lymphoblastoid cell lines were obtained. 422 of this samples were sequenced within the 1000 Genome Project Phase 1. Here, we make available a subset of bi-allelic SNPs and transcript expected counts for CEPH (CEU) population that corresponds to 50 randomly selected genes from chromosome 19.
#' 
#' @format 
#' \code{dataSQTL_counts} contains transcript expected counts. A data.frame with 93 variables:
#' \itemize{
#'   \item \code{gene_id}: Gene IDs.
#'   \item \code{transcript_id}: Transcript IDs.
#'   \item \code{NA06984, NA06985, ..., NA12890}: Quantification of transcripts for 91 individuals.
#'   }
#'   
#'   
#'  \code{dataSQTL_gene_ranges} with the information about gene location. A \code{\linkS4class{GRanges}} object with gene ranges.
#'  
#'  
#' \code{dataSQTL_genotypes} contains genotypes coded as follows: 0 for ref/ref, 1 for ref/not ref, 2 for not ref/not ref, -1 or NA for missing value. A data.frame with 95 variables:
#' \itemize{
#'   \item \code{chr}: Chromosome IDs.
#'   \item \code{start}: Start of SNP location.
#'   \item \code{end}: End of SNP location.
#'   \item \code{snp_id}: SNP IDs.
#'   \item \code{NA06984, NA06985, ..., NA12890}: Genotypes for 91 individuals.
#' }
#' 
#' 
#' \code{dataSQTL_dmSQTLdata} is a \code{DM} package object created through the sQTL analysis pipeline. 
#' 
#' @source 
#' Lappalainen T, Sammeth M, Friedländer MR, et al. Transcriptome and genome sequencing uncovers functional variation in humans. Nature. 2013;501(7468):506–11.
#' 
#' Genotypes and transcript quantification were downloaded from http://www.ebi.ac.uk/Tools/geuvadis-das/.
#' 
#' Gene annotation Gencode v12 from http://www.gencodegenes.org/releases/12.html.
"dataSQTL_counts"



#' @rdname dataSQTL_counts
"dataSQTL_genotypes"

#' @rdname dataSQTL_counts
"dataSQTL_gene_ranges"

#' @rdname dataSQTL_counts
"dataSQTL_dmSQTLdata"











