#' @import BiocGenerics
#' @import methods 
#' @import BiocParallel
#' @import edgeR
#' @import GenomicRanges
#' @import S4Vectors
#' @import ggplot2
#' @importFrom reshape2 melt
#' @import VennDiagram
#' @import grid
NULL





################################################################################
### Data documentation
################################################################################



#' Sample data for differential splicing analysis
#'
#' A subset (100 randomly selected genes) of exonic bin counts computed with DEXSeq package in a simulation study based on Drosophila genome. In this study, data was generated to mimic an RNA-seq assay for 6 samples separated into 2 conditions (3 versus 3 samples). Differential transcript usage was induced for 1000 genes by swapping, between the 2 conditions, the expression of two most abundant transcripts. Our subset consists of 20 genes with differential splicing (status 1) and 80 without (status 0). For more details about this simulations, see the reference in Source paragraph.
#' 
#' @format 
#' \code{dataDS_counts} contains exonic bin counts. A data frame with 7 variables:
#' \itemize{
#'   \item \code{group_id}: Exonic bin IDs, which contain gene IDs and bin numbers separated with ':'.
#'   \item \code{sample_1, sample_2, ...}: Quantification of the exonic bins in 6 samples.
#'   }
#'   
#' \code{dataDS_metadata} describes the samples. A data frame with 2 variables:
#' \itemize{
#'   \item \code{sample_id}: Sample IDs.
#'   \item \code{group}: Grouping into two conditions.
#' }
#' 
#' \code{dataDS_status} contains the information about the differential splicing status of genes. A data frame with 2 variables:
#' \itemize{
#'   \item \code{gene_id}: Gene IDs.
#'   \item \code{status}: 1 for genes with differential splicing and 0 for genes without DS.
#' }
#' 
#' \code{dataDS_dmDSdata}, \code{dataDS_dmDSdispersion}, \code{dataDS_dmDSfit} and \code{dataDS_dmDStest} are \code{DM} package objects created through the differential splicing analysis pipeline. See Examples.
#' 
#' @source 
#' Soneson C, Matthes KL, Nowicka M, Law CW, Robinson MD. Differential transcript usage from RNA-seq data : isoform pre-filtering improves performance of count-based methods. 2015.
#' 
#' \url{http://imlspenticton.uzh.ch/robinson_lab/splicing_comparison/}
#' 
#' @examples 
#' 
#' ### Differential splicing analysis
#' \donttest{
#' 
#'  ### Create dmDSdata object
#' 
#'  # counts
#'  head(dataDS_counts)
#'  # metadata
#'  head(dataDS_metadata)
#'  
#'  group_split <- limma::strsplit2(dataDS_counts[, 1], ":")
#'  
#'  d <- dmDSdata(counts = dataDS_counts[, -1], gene_id = group_split[, 1], 
#'    feature_id = group_split[, 2], sample_id = dataDS_metadata$sample_id, 
#'    group = dataDS_metadata$group)
#'  
#'  plotData(d)
#'  
#'  ## dataDS_dmDSdata <- d
#'  
#'  ### Filtering
#'  
#'  # Check what is the minimal number of replicates per condition 
#'  table(samples(d)$group)
#'  
#'  d <- dmFilter(d, min_samps_gene_expr = 3, min_samps_feature_prop = 3)
#'  plotData(d)
#'  
#'  ### Calculate dispersion
#'  
#'  # If possible, increase the number of workers
#'  d <- dmDispersion(d, BPPARAM = BiocParallel::MulticoreParam(workers = 1))
#' 
#'  plotDispersion(d)
#'  
#'  ## dataDS_dmDSdispersion <- d
#'  
#'  ### Fit full model proportions
#'  
#'  # If possible, increase the number of workers
#'  d <- dmFit(d, BPPARAM = BiocParallel::MulticoreParam(workers = 1))
#'  
#'  ## dataDS_dmDSfit <- d
#'  
#'  ### Fit null model proportions and test for DS
#'  
#'  # If possible, increase the number of workers
#'  d <- dmTest(d, BPPARAM = BiocParallel::MulticoreParam(workers = 1))
#' 
#'  plotTest(d)
#' 
#'  ## dataDS_dmDStest <- d
#' 
#'  ### Plot feature proportions for top DS gene
#' 
#'  res <- results(d)
#'  res <- res[order(res$pvalue, decreasing = FALSE), ]
#' 
#'  gene_id <- res$gene_id[1]
#' 
#'  plotFit(d, gene_id = gene_id)
#'  plotFit(d, gene_id = gene_id, plot_type = "lineplot")
#'  plotFit(d, gene_id = gene_id, plot_type = "ribbonplot")
#'
#' }
#' 
#' 
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
"dataDS_dmDSfit"

#' @rdname dataDS_counts
"dataDS_dmDStest"







#' Sample data for sQTL analysis
#'
#' A subset of data from GEUVADIS project where 465 RNA-seq samples from lymphoblastoid cell lines were obtained. 422 of this samples were sequenced within the 1000 Genome Project Phase 1. Here, we make available a subset of bi-allelic SNPs and transcript expected counts for CEPH (CEU) population that corresponds to 50 randomly selected genes from chromosome 19.
#' 
#' @format 
#' \code{dataSQTL_counts} contains transcript expected counts. A data frame with 93 variables:
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
#' \code{dataSQTL_genotypes} contains genotypes coded as follows: 0 for ref/ref, 1 for ref/not ref, 2 for not ref/not ref, -1 or NA for missing value. A data frame with 95 variables:
#' \itemize{
#'   \item \code{chr}: Chromosome IDs.
#'   \item \code{start}: Start of SNP location.
#'   \item \code{end}: End of SNP location.
#'   \item \code{snp_id}: SNP IDs.
#'   \item \code{NA06984, NA06985, ..., NA12890}: Genotypes for 91 individuals.
#' }
#' 
#' 
#' \code{dataSQTL_dmSQTLdata} and \code{dataSQTL_dmSQTLtest} are \code{DM} package objects created through the sQTL analysis pipeline. See Examples.
#' 
#' @source 
#' Lappalainen T, Sammeth M, Friedlander MR, et al. Transcriptome and genome sequencing uncovers functional variation in humans. Nature. 2013;501(7468):506-11.
#' 
#' Genotypes and transcript quantification were downloaded from http://www.ebi.ac.uk/Tools/geuvadis-das/.
#' 
#' Gene annotation Gencode v12 from http://www.gencodegenes.org/releases/12.html.
#' 
#' @examples 
#' 
#' ### sQTL analysis
#' \donttest{
#' 
#'  ### Create dmSQTLdata object
#'  
#'  # counts
#'  head(dataSQTL_counts)
#'  # gene_ranges
#'  dataSQTL_gene_ranges
#'  # genotypes 
#'  head(dataSQTL_genotypes)
#'  
#'  ## gene_ranges with names!
#'  gene_ranges <- dataSQTL_gene_ranges
#'  names(gene_ranges) <- S4Vectors::mcols(gene_ranges)$name
#'  
#'  ## snp_ranges with names!
#'  snp_ranges <- GenomicRanges::GRanges(S4Vectors::Rle(dataSQTL_genotypes$chr), 
#'    IRanges::IRanges(dataSQTL_genotypes$start, dataSQTL_genotypes$end))
#'  names(snp_ranges) <- dataSQTL_genotypes$snp_id 
#'  
#'  ## Check if samples in count and genotypes are in the same order
#'  all(colnames(dataSQTL_counts[, -(1:2)]) == colnames(dataSQTL_genotypes[, -(1:4)]))
#'  sample_id <- colnames(dataSQTL_counts[, -(1:2)])
#'  
#'  d <- dmSQTLdataFromRanges(counts = dataSQTL_counts[, -(1:2)], 
#'    gene_id = dataSQTL_counts$gene_id, feature_id = dataSQTL_counts$transcript_id, 
#'    gene_ranges = gene_ranges, genotypes = dataSQTL_genotypes[, -(1:4)], 
#'    snp_id = dataSQTL_genotypes$snp_id, snp_ranges = snp_ranges, 
#'    sample_id = sample_id, window = 5e3, 
#'    BPPARAM = BiocParallel::MulticoreParam(workers = 1))
#'  
#'  plotData(d)
#'  
#'  ## dataSQTL_dmSQTLdata <- d
#'  
#'  ### Filtering
#'  
#'  # If possible, increase the number of workers
#'  d <- dmFilter(d, min_samps_gene_expr = 70, min_samps_feature_prop = 5,
#'    minor_allele_freq = 5, BPPARAM = BiocParallel::MulticoreParam(workers = 1))
#'  plotData(d)
#'  
#'  ### Calculate dispersion
#'  
#'  # If possible, increase the number of workers
#'  d <- dmDispersion(d, BPPARAM = BiocParallel::MulticoreParam(workers = 1))
#' 
#'  plotDispersion(d)
#'  
#'  ### Fit full model proportions
#'  
#'  # If possible, increase the number of workers
#'  d <- dmFit(d, BPPARAM = BiocParallel::MulticoreParam(workers = 1))
#' 
#'  ### Fit null model proportions and test for sQTLs
#'  
#'  # If possible, increase the number of workers
#'  d <- dmTest(d, BPPARAM = BiocParallel::MulticoreParam(workers = 1))
#' 
#'  plotTest(d)
#' 
#'  ## dataSQTL_dmSQTLtest <- d
#' 
#'  ### Plot feature proportions for top sQTL
#'  
#'  res <- results(d)
#'  res <- res[order(res$pvalue, decreasing = FALSE), ]
#' 
#'  gene_id <- res$gene_id[1]
#'  snp_id <- res$snp_id[1]
#' 
#'  plotFit(d, gene_id, snp_id)
#'  plotFit(d, gene_id, snp_id, plot_type = "boxplot2", order = FALSE)
#'  plotFit(d, gene_id, snp_id, plot_type = "ribbonplot")
#'  
#' }
#' 
"dataSQTL_counts"



#' @rdname dataSQTL_counts
"dataSQTL_genotypes"

#' @rdname dataSQTL_counts
"dataSQTL_gene_ranges"

#' @rdname dataSQTL_counts
"dataSQTL_dmSQTLdata"

#' @rdname dataSQTL_counts
"dataSQTL_dmSQTLtest"








