# Prepare data for examples and vignette 


# library(devtools); library(GenomicRanges); library(BiocParallel); library(edgeR)


# Rfiles <- list.files("R/", full.names=TRUE); for(i in Rfiles) source(i)

library(DM)


## Input files: transcript expression, gene location and genotype information
data_dir <- "data-raw/geuvadis/data/"


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

# library(rtracklayer)

genes_path = paste0(data_dir, "annotation/genes_noChr.bed")
gene_ranges = rtracklayer::import(genes_path)
names(gene_ranges) <- S4Vectors::mcols(gene_ranges)$name


sample_id <- colnames(genotypes)

window <- 5e3




data_org <- dmSQTLdataFromRanges(counts, gene_id_counts, feature_id_counts, gene_ranges, genotypes, snp_id_genotypes, snp_ranges, sample_id, window = 5e3)


set.seed(1)
genes_subset <- names(data_org@counts)[sample(length(names(data_org@counts)), size = 1000, replace = FALSE)]



data_sub <- new("dmSQTLdata", counts = data_org@counts[genes_subset], genotypes = data_org@genotypes[genes_subset], samples = data_org@samples)

dmSQTLplotData(data_sub)



data_sub_filt_org <- data_sub_filt <- dmSQTLfilter(data_sub, min_samps_gene_expr = 70, min_gene_expr = 1, min_samps_feature_prop = 5, min_feature_prop = 0.1, max_features = Inf, minor_allel_freq = 0.05, BPPARAM = BiocParallel::MulticoreParam(workers = 3))

dmSQTLplotData(data_sub_filt)


### This is the data to be used in examples 

genes_subset <- names(data_sub_filt@genotypes@partitioning)[order(IRanges::width(data_sub_filt@genotypes@partitioning), decreasing = FALSE)[1:50]]

data_sub <- new("dmSQTLdata", counts = data_org@counts[genes_subset], genotypes = data_org@genotypes[genes_subset], samples = data_org@samples)

dmSQTLplotData(data_sub)


### Get the basic components 


## counts
gene_id_counts <- rep(names(data_sub@counts), width(data_sub@counts@partitioning))
feature_id_counts <- rownames(data_sub@counts@unlistData)
counts <- data_sub@counts@unlistData

dataSQTL_counts <- data.frame(group_id = paste0(gene_id_counts, ":", feature_id_counts), counts, row.names = NULL, stringsAsFactors = FALSE)

use_data(dataSQTL_counts, overwrite = TRUE)


## gene_ranges
dataSQTL_gene_ranges <- gene_ranges[names(data_sub@counts)]
names(dataSQTL_gene_ranges) <- NULL

use_data(dataSQTL_gene_ranges, overwrite = TRUE)


## genotypes

snp_id_genotypes <- unique(rownames(data_sub@genotypes@unlistData))

dataSQTL_genotypes <- data.frame(genotypes_raw[genotypes_raw$snpId %in% snp_id_genotypes, ], row.names = NULL, stringsAsFactors = FALSE)

colnames(dataSQTL_genotypes)[4] <- "snp_id"

use_data(dataSQTL_genotypes, overwrite = TRUE)


### Use the basic components to create dmSQTLdata object

### counts
head(dataSQTL_counts)
counts <- as.matrix(dataSQTL_counts[, -1])

group_id <- dataSQTL_counts[,1]
group_split <- limma::strsplit2(group_id, ":")
gene_id_counts <- group_split[, 1]
feature_id_counts <- group_split[, 2]

### gene_ranges
dataSQTL_gene_ranges
gene_ranges <- dataSQTL_gene_ranges
names(gene_ranges) <- S4Vectors::mcols(gene_ranges)$name


### genotypes
head(dataSQTL_genotypes)
genotypes <- as.matrix(dataSQTL_genotypes[, -(1:4)])

snp_id_genotypes <- dataSQTL_genotypes$snp_id

snp_ranges <- GenomicRanges::GRanges(S4Vectors::Rle(dataSQTL_genotypes$chr), IRanges::IRanges(dataSQTL_genotypes$start, dataSQTL_genotypes$end))
names(snp_ranges) <- dataSQTL_genotypes$snp_id

all(colnames(counts) == colnames(genotypes))

sample_id <- colnames(counts)

### create dmSQTLdata object 
data <- dmSQTLdataFromRanges(counts, gene_id_counts, feature_id_counts, gene_ranges, genotypes, snp_id_genotypes, snp_ranges, sample_id, window = 5e3)

dmSQTLplotData(data)

dataSQTL_dmSQTLdata <- data

use_data(dataSQTL_dmSQTLdata, overwrite = TRUE)



### filtering

data <- dataSQTL_dmSQTLdata
dmSQTLplotData(data)
data <- dmSQTLfilter(data)
dmSQTLplotData(data)


### dispersion

data <- dataSQTL_dmSQTLdata
dmSQTLplotData(data)

data <- dmSQTLfilter(data)
dmSQTLplotData(data)


data <- dmSQTLdispersion(data)

dataSQTL_dmSQTLdispersion <- data

use_data(dataSQTL_dmSQTLdispersion, overwrite = TRUE)


dmSQTLplotDispersion(data)


### fit

data <- dataSQTL_dmSQTLdispersion

data <- dmSQTLfit(data)


snp_id <- "snp_19_12796435"
gene_id <- "ENSG00000132004.7"

dmSQTLplotFit(data, gene_id, snp_id)


### test

data <- dmSQTLtest(data)

dmSQTLplotTest(data)


snp_id <- "snp_19_48981946"
gene_id <- "ENSG00000105443.8"

dmSQTLplotFit(data, gene_id, snp_id)










