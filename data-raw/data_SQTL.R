# Prepare data for examples and vignette 

setwd("/home/gosia/R/multinomial_project/package_devel/DM/data-raw/geuvadis/")

library(DM)

library(devtools)


########################################################
# load raw data
########################################################

## Input files: transcript expression, gene location and genotype information
data_dir <- "data/"


### read counts 
counts_path <- paste0(data_dir, "expression/trExpCount_CEU.tsv")
counts_raw <- read.table(counts_path, header = TRUE, as.is = TRUE)

gene_id <- counts_raw$geneId
feature_id <- counts_raw$trId
counts <- as.matrix(counts_raw[, -c(1:2)])


### read genotypes
chr = "19"
genotypes_raw <- read.table(paste0(data_dir, "genotypes/snps_CEU_chr", chr ,".tsv"), header = TRUE, sep = "\t", as.is = TRUE)

snp_ranges <- GenomicRanges::GRanges(S4Vectors::Rle(chr, nrow(genotypes_raw)), IRanges::IRanges(genotypes_raw$start, genotypes_raw$end))
names(snp_ranges) <- genotypes_raw$snpId
snp_id <- genotypes_raw$snpId

genotypes <- as.matrix(genotypes_raw[, -c(1:4)])


### read ranges
genes_path = paste0(data_dir, "annotation/genes_noChr.bed")
gene_ranges = rtracklayer::import(genes_path)
names(gene_ranges) <- S4Vectors::mcols(gene_ranges)$name


sample_id <- colnames(genotypes)

window <- 5e3


########################################################
# Find subset
########################################################


do <- dmSQTLdataFromRanges(counts, gene_id, feature_id, gene_ranges, genotypes, snp_id, snp_ranges, sample_id, window = 5e3)

plotData(do, out_dir = "./")

dof <- dmFilter(do, min_samps_gene_expr = 70, min_gene_expr = 1, min_samps_feature_prop = 5, min_feature_prop = 0.1, max_features = Inf, minor_allel_freq = 5, BPPARAM = BiocParallel::MulticoreParam(workers = 10))

plotData(dof, out_dir = "./")

oo <- order(width(dof@genotypes), decreasing = FALSE)

genes_subset <- names(dof)[oo][1:50]


ds <- do[genes_subset, ]

plotData(ds, out_dir = "./")

########################################################
# Get the basic components
########################################################

## counts
gene_id <- rep(names(ds@counts), width(ds@counts))
feature_id <- rownames(ds@counts)
counts <- ds@counts@unlistData

dataSQTL_counts <- data.frame(gene_id = gene_id, transcript_id = feature_id, counts, row.names = NULL, stringsAsFactors = FALSE)

use_data(dataSQTL_counts, overwrite = TRUE)


## gene_ranges
dataSQTL_gene_ranges <- gene_ranges[genes_subset]
names(dataSQTL_gene_ranges) <- NULL

use_data(dataSQTL_gene_ranges, overwrite = TRUE)


## genotypes

snp_id <- unique(rownames(ds@genotypes))

dataSQTL_genotypes <- data.frame(genotypes_raw[genotypes_raw$snpId %in% snp_id, ], row.names = NULL, stringsAsFactors = FALSE)
colnames(dataSQTL_genotypes)[4] <- "snp_id"


use_data(dataSQTL_genotypes, overwrite = TRUE)

########################################################
### Use the basic components to create dmSQTLdata object
########################################################

### counts
head(dataSQTL_counts)

counts <- dataSQTL_counts[, -(1:2)]
gene_id <- group_split[, 1]
feature_id <- group_split[, 2]

### gene_ranges
dataSQTL_gene_ranges
gene_ranges <- dataSQTL_gene_ranges
names(gene_ranges) <- S4Vectors::mcols(gene_ranges)$name

### genotypes
head(dataSQTL_genotypes)
genotypes <- dataSQTL_genotypes[, -(1:4)]

snp_id <- dataSQTL_genotypes$snp_id

snp_ranges <- GenomicRanges::GRanges(S4Vectors::Rle(dataSQTL_genotypes$chr), IRanges::IRanges(dataSQTL_genotypes$start, dataSQTL_genotypes$end))
names(snp_ranges) <- dataSQTL_genotypes$snp_id

all(colnames(counts) == colnames(genotypes))

sample_id <- colnames(counts)

### create dmSQTLdata object 
d <- dmSQTLdataFromRanges(counts, gene_id, feature_id, gene_ranges, genotypes, snp_id, snp_ranges, sample_id, window = 5e3)


dataSQTL_dmSQTLdata <- d

use_data(dataSQTL_dmSQTLdata, overwrite = TRUE)


d <- dataSQTL_dmSQTLdata
plotData(d)

dev.off()




### filtering

dd <- dataSQTL_dmSQTLdata
plotData(dd)
d <- dmSQTLfilter(d)
plotData(d)




### dispersion

data <- dataSQTL_dmSQTLdata
plotData(data)

data <- dmFilter(data)
plotData(data)


data <- dmDispersion(data)

dataSQTL_dmSQTLdispersion <- data

use_data(dataSQTL_dmSQTLdispersion, overwrite = TRUE)


plotDispersion(data)


### fit

data <- dataSQTL_dmSQTLdispersion

data <- dmFit(data)


snp_id <- "snp_19_12796435"
gene_id <- "ENSG00000132004.7"

plotFit(data, gene_id, snp_id)


### test

data <- dmLRT(data)

plotLRT(data)


snp_id <- "snp_19_48981946"
gene_id <- "ENSG00000105443.8"

plotFit(data, gene_id, snp_id)










