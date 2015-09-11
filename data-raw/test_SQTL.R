# Test DM on entire data set (chr19)

library(DM)

setwd("/home/gosia/R/multinomial_project/package_devel/DM/data-raw/geuvadis/")


########################################################
# sqtl data
########################################################

## Input files: transcript expression, gene location and genotype information
data_dir <- "data/"


### read counts 
counts_path <- paste0(data_dir, "expression/trExpCount_CEU.tsv")
counts_raw <- read.table(counts_path, header = TRUE, as.is = TRUE)

gene_id_counts <- counts_raw$geneId
feature_id_counts <- counts_raw$trId
counts <- as.matrix(counts_raw[, -c(1:2)])


### read genotypes
chr = "19"
genotypes_raw <- read.table(paste0(data_dir, "genotypes/snps_CEU_chr", chr ,".tsv"), header = TRUE, sep = "\t", as.is = TRUE)

snp_ranges <- GenomicRanges::GRanges(S4Vectors::Rle(chr, nrow(genotypes_raw)), IRanges::IRanges(genotypes_raw$start, genotypes_raw$end))
names(snp_ranges) <- genotypes_raw$snpId
snp_id_genotypes <- genotypes_raw$snpId

genotypes <- as.matrix(genotypes_raw[, -c(1:4)])


### read ranges

genes_path = paste0(data_dir, "annotation/genes_noChr.bed")
gene_ranges = rtracklayer::import(genes_path)
names(gene_ranges) <- S4Vectors::mcols(gene_ranges)$name


sample_id <- colnames(genotypes)

window <- 5e3


d <- dmSQTLdataFromRanges(counts, gene_id_counts, feature_id_counts, gene_ranges, genotypes, snp_id_genotypes, snp_ranges, sample_id, window = 5e3)


all(names(d@counts) == names(d@genotypes))


d <- dmFilter(d, BPPARAM = BiocParallel::MulticoreParam(workers = 10))


plotData(d, "./")


# save(d, file = "data-raw/d.Rdata")



ds <- dmDispersion(d, BPPARAM = BiocParallel::MulticoreParam(workers = 10))


# save(ds, file = "data-raw/ds.Rdata")


plotDispersion(ds, "./")


df <- dmFit(ds, BPPARAM = BiocParallel::MulticoreParam(workers = 10))


snp_id <- "snp_19_54704760"
gene_id <- "ENSG00000170889.9"

plotFit(df, gene_id, snp_id, out_dir = "./")



dt <- dmLRT(df, BPPARAM = BiocParallel::MulticoreParam(workers = 10))


plotLRT(dt, out_dir = "./")


plotFit(dt, gene_id, snp_id, out_dir = "./", plot_type = "boxplot2")






