#' @import BiocGenerics
#' @import methods 
#' @import BiocParallel
#' @import edgeR
#' @import GenomicRanges
#' @import S4Vectors
#' @import ggplot2
#' @importFrom reshape2 melt
#' @import grid
NULL





################################################################################
### Data documentation
################################################################################



#' Sample data for differential splicing analysis
#'
#' We use a subset of HTSeq exonic bin counts from \code{pasilla} package.
#' 
#' @format 
#' \code{data_dmDSdata} is a \code{\linkS4class{dmDSdata}} object. See Examples.
#' 
#' @source 
#' Brooks AN, Yang L, Duff MO, et al. Conservation of an RNA regulatory map between Drosophila and mammals. Genome Res. 2011;21(2):193-202.
#' 
#' \code{pasilla} package.
#' 
#' @return 
#' \code{data_dmDSdata}
#' 
#' @examples 
#' 
#' #############################
#' ### Create dmDSdata object
#' #############################
#' ### Get HTSeq exonic bin counts from 'pasilla' package
#' 
#' library(pasilla)
#' 
#' data_dir  <- system.file("extdata", package="pasilla")
#' count_files <- list.files(data_dir, pattern="fb.txt$", full.names=TRUE)
#' count_files
#' 
#' # Create a data frame with htseq counts
#' htseq_list <- lapply(1:length(count_files), function(i){
#'   # i = 1
#'   htseq <- read.table(count_files[i], header = FALSE, as.is = TRUE)
#'   colnames(htseq) <- c("group_id", gsub("fb.txt", "", strsplit(count_files[i], 
#'      "extdata/")[[1]][2]))
#'   return(htseq)
#' })
#' 
#' htseq_counts <- Reduce(function(...) merge(..., by = "group_id", all=TRUE, 
#'    sort = FALSE), htseq_list)
#' tail(htseq_counts)
#' htseq_counts <- htseq_counts[!grepl(pattern = "_", htseq_counts$group_id), ]
#' 
#' group_split <- limma::strsplit2(htseq_counts[, 1], ":")
#' 
#' d <- dmDSdata(counts = htseq_counts[, -1], gene_id = group_split[, 1], 
#'    feature_id = group_split[, 2], sample_id = colnames(htseq_counts)[-1], 
#'    group = gsub("[1-4]", "", colnames(htseq_counts)[-1]))
#' 
#' plotData(d)
#' 
#' # Use a subset of genes, which is defined in the following file
#' genes_subset = readLines(file.path(data_dir, "geneIDsinsubset.txt"))
#' d <- d[names(d) %in% genes_subset, ]
#' 
#' plotData(d)
#' 
#' ## data_dmDSdata <- d
#' 
#' ###################################
#' ### Differential splicing analysis
#' ###################################
#' # If possible, increase the number of workers in BPPARAM
#' 
#' d <- data_dmDSdata
#' 
#' head(counts(d))
#' samples(d)
#' head(names(d))
#' length(d)
#' d[1:20, ]
#' d[1:20, 1:3]
#' 
#' ### Filtering
#' # Check what is the minimal number of replicates per condition 
#' table(samples(d)$group)
#' d <- dmFilter(d, min_samps_gene_expr = 6, min_samps_feature_expr = 3, 
#'  min_samps_feature_prop = 3)
#' plotData(d)
#' 
#' ### Calculate dispersion
#' d <- dmDispersion(d, BPPARAM = BiocParallel::MulticoreParam(workers = 3))
#' plotDispersion(d)
#' 
#' head(mean_expression(d))
#' common_dispersion(d)
#' head(genewise_dispersion(d))
#' 
#' ### Fit full model proportions
#' d <- dmFit(d, BPPARAM = BiocParallel::MulticoreParam(workers = 1))
#' 
#' head(proportions(d))
#' head(statistics(d))
#' 
#' ### Fit null model proportions and test for DS
#' d <- dmTest(d, BPPARAM = BiocParallel::MulticoreParam(workers = 1))
#' plotTest(d)
#' 
#' head(proportions(d))
#' head(statistics(d))
#' head(results(d))
#' 
#' ### Plot feature proportions for top DS gene
#' res <- results(d)
#' res <- res[order(res$pvalue, decreasing = FALSE), ]
#' 
#' gene_id <- res$gene_id[1]
#' 
#' plotFit(d, gene_id = gene_id)
#' plotFit(d, gene_id = gene_id, plot_type = "lineplot")
#' plotFit(d, gene_id = gene_id, plot_type = "ribbonplot")
#' 
#' 
"data_dmDSdata"






#' Sample data for sQTL analysis
#'
#' A subset of data from GEUVADIS project where 462 RNA-Seq samples from lymphoblastoid cell lines were obtained. The genome sequencing data of the same individuals is provided by the 1000 Genomes Project. The samples in this project come from five populations: CEPH (CEU), Finns (FIN), British (GBR), Toscani (TSI) and Yoruba (YRI). Here, we make available subsets of bi-allelic SNPs and transcript expected counts for CEPH population (91 individuals) that correspond to 50 randomly selected genes from chromosome 19. For the details on how this data was preprocessed, see the vignette.
#' 
#' @format 
#' \code{data_dmSQTLdata} is a \code{\linkS4class{dmSQTLdata}} object. See Examples.
#' 
#' @source 
#' Lappalainen T, Sammeth M, Friedlander MR, et al. Transcriptome and genome sequencing uncovers functional variation in humans. Nature. 2013;501(7468):506-11.
#' 
#' Genotypes and transcript quantification were downloaded from http://www.ebi.ac.uk/Tools/geuvadis-das/.
#' 
#' Gene annotation Gencode v12 from http://www.gencodegenes.org/releases/12.html.
#' 
#' @return 
#' \code{data_dmSQTLdata}
#' 
#' @examples 
#' 
#' #############################
#' ### Create dmSQTLdata object
#' #############################
#' 
#' library(GenomicRanges)
#' library(rtracklayer)
#' 
#' data_dir  <- system.file("extdata", package = "DRIMSeq")
#' 
#' 
#' # gene_ranges with names!
#' gene_ranges <- import(paste0(data_dir, "/genes_subset.bed"))
#' names(gene_ranges) <- mcols(gene_ranges)$name
#' 
#' counts <- read.table(paste0(data_dir, "/TrQuantCount_CEU_subset.tsv"), 
#'    header = TRUE, sep = "\t", as.is = TRUE)
#' 
#' genotypes <- read.table(paste0(data_dir, "/genotypes_CEU_subset.tsv"), 
#'    header = TRUE, sep = "\t", as.is = TRUE)
#' 
#' # snp_ranges with names!
#' snp_ranges <- GRanges(Rle(genotypes$chr), IRanges(genotypes$start, 
#'    genotypes$end))
#' names(snp_ranges) <- genotypes$snpId 
#' 
#' ## Check if samples in count and genotypes are in the same order
#' all(colnames(counts[, -(1:2)]) == colnames(genotypes[, -(1:4)]))
#' sample_id <- colnames(counts[, -(1:2)])
#' 
#' 
#' d <- dmSQTLdataFromRanges(counts = counts[, -(1:2)], gene_id = counts$geneId, 
#'    feature_id = counts$trId, gene_ranges = gene_ranges, 
#'    genotypes = genotypes[, -(1:4)], snp_id = genotypes$snpId, 
#'    snp_ranges = snp_ranges, sample_id = sample_id, window = 5e3, 
#'    BPPARAM = BiocParallel::MulticoreParam(workers = 1))
#' 
#' plotData(d)
#' 
#' 
#' ## data_dmSQTLdata <- d
#' 
#' #############################
#' ### sQTL analysis
#' #############################
#' # If possible, increase the number of workers in BPPARAM
#' 
#' d <- data_dmSQTLdata
#' 
#' head(names(d))
#' length(d)
#' d[1:10, ]
#' d[1:10, 1:10]
#' 
#' ### Filtering
#' d <- dmFilter(d, min_samps_gene_expr = 70, min_samps_feature_expr = 5, 
#'    min_samps_feature_prop = 5, minor_allele_freq = 5, 
#'    BPPARAM = BiocParallel::MulticoreParam(workers = 1))
#' plotData(d)
#' 
#' 
#' ### Calculate dispersion
#' d <- dmDispersion(d, BPPARAM = BiocParallel::MulticoreParam(workers = 3))
#' plotDispersion(d)
#' 
#' 
#' ### Fit full model proportions
#' d <- dmFit(d, BPPARAM = BiocParallel::MulticoreParam(workers = 1))
#' 
#' 
#' ### Fit null model proportions and test for sQTLs
#' d <- dmTest(d, BPPARAM = BiocParallel::MulticoreParam(workers = 1))
#' plotTest(d)
#' 
#' head(results(d))
#' 
#' ### Plot feature proportions for top sQTL
#' res <- results(d)
#' res <- res[order(res$pvalue, decreasing = FALSE), ]
#' 
#' gene_id <- res$gene_id[1]
#' snp_id <- res$snp_id[1]
#' 
#' plotFit(d, gene_id, snp_id)
#' plotFit(d, gene_id, snp_id, plot_type = "boxplot2", order = FALSE)
#' plotFit(d, gene_id, snp_id, plot_type = "ribbonplot")
#' 
"data_dmSQTLdata"











