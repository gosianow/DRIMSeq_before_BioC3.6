# Prepare data for examples and vignette 

# setwd("/home/gosia/R/multinomial_project/package_devel/DM/data-raw/geuvadis/")

library(DM)

library(devtools)


########################################################
# load raw data
########################################################

## Input files: transcript expression, gene location and genotype information
data_dir <- "data-raw/geuvadis/data/"


### read counts 
counts_path <- paste0(data_dir, "expression/trExpCount_CEU.tsv")

counts_raw <- readLines(counts_path)
tools::showNonASCII(counts_raw)


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
# Find subset - genes with little number of snps
########################################################


do <- dmSQTLdataFromRanges(counts, gene_id, feature_id, gene_ranges, genotypes, snp_id, snp_ranges, sample_id, window = 5e3, BPPARAM = BiocParallel::MulticoreParam(workers = 2))

plotData(do, out_dir = "./")

dof <- dmFilter(do, min_samps_gene_expr = 70, min_gene_expr = 1, min_samps_feature_prop = 5, min_feature_prop = 0.1, max_features = Inf, minor_allele_freq = 5, BPPARAM = BiocParallel::MulticoreParam(workers = 2))

plotData(dof, out_dir = "./")


### take 50 genes
oo <- order(width(dof@genotypes), decreasing = FALSE)

genes_subset <- names(dof)[oo][1:50]


ds <- do[genes_subset, ]

plotData(ds, out_dir = "./")






########################################################
# Get the basic components
########################################################

## counts
gene_id <- as.character(rep(names(ds@counts), width(ds@counts)))
feature_id <- as.character(rownames(ds@counts))
counts <- ds@counts@unlistData

tools::showNonASCII(gene_id)
tools::showNonASCII(feature_id)

# gi <- iconv(gene_id, from = "ASCII", to = "ASCII", sub = "")

gene_id <- gsub("[^\x20-\x7E]", "", gene_id)
feature_id <- gsub("[^\x20-\x7E]", "", feature_id)


dataSQTL_counts <- data.frame(gene_id = gene_id, transcript_id = feature_id, counts, row.names = NULL, stringsAsFactors = FALSE)

tools::showNonASCII(colnames(dataSQTL_counts))

for(i in 1:ncol(dataSQTL_counts))
  print(tools::showNonASCII(dataSQTL_counts[, i]))


use_data(dataSQTL_counts, overwrite = TRUE)


## gene_ranges
dataSQTL_gene_ranges <- gene_ranges[genes_subset]
names(dataSQTL_gene_ranges) <- NULL

use_data(dataSQTL_gene_ranges, overwrite = TRUE)


## genotypes

snp_id <- unique(ds@blocks@unlistData[, "snp_id"])

dataSQTL_genotypes <- data.frame(genotypes_raw[genotypes_raw$snpId %in% snp_id, ], row.names = NULL, stringsAsFactors = FALSE)
colnames(dataSQTL_genotypes)[4] <- "snp_id"


use_data(dataSQTL_genotypes, overwrite = TRUE)

########################################################
### Use the basic components to create dmSQTLdata object
########################################################



### Create dmSQTLdata object

# counts
head(dataSQTL_counts)
# gene_ranges
dataSQTL_gene_ranges
# genotypes 
head(dataSQTL_genotypes)

# gene_ranges with names!
gene_ranges <- dataSQTL_gene_ranges
names(gene_ranges) <- S4Vectors::mcols(gene_ranges)$name

# snp_ranges with names!
snp_ranges <- GenomicRanges::GRanges(S4Vectors::Rle(dataSQTL_genotypes$chr), 
                                     IRanges::IRanges(dataSQTL_genotypes$start, dataSQTL_genotypes$end))
names(snp_ranges) <- dataSQTL_genotypes$snp_id 

## Check if samples in count and genotypes are in the same order
all(colnames(dataSQTL_counts[, -(1:2)]) == colnames(dataSQTL_genotypes[, -(1:4)]))
sample_id <- colnames(dataSQTL_counts[, -(1:2)])

d <- dmSQTLdataFromRanges(counts = dataSQTL_counts[, -(1:2)], 
                          gene_id = dataSQTL_counts$gene_id, feature_id = dataSQTL_counts$transcript_id, 
                          gene_ranges = gene_ranges, genotypes = dataSQTL_genotypes[, -(1:4)], 
                          snp_id = dataSQTL_genotypes$snp_id, snp_ranges = snp_ranges, 
                          sample_id = sample_id, window = 5e3, 
                          BPPARAM = BiocParallel::MulticoreParam(workers = 1))

plotData(d)

dataSQTL_dmSQTLdata <- d



use_data(dataSQTL_dmSQTLdata, overwrite = TRUE)



### Filtering

d <- dataSQTL_dmSQTLdata

# If possible, increase the number of workers
d <- dmFilter(d, BPPARAM = BiocParallel::MulticoreParam(workers = 2))

plotData(d)


### Calculate dispersion

# If possible, increase the number of workers
d <- dmDispersion(d, BPPARAM = BiocParallel::MulticoreParam(workers = 2))

plotDispersion(d)

### Fit full model proportions

# If possible, increase the number of workers
d <- dmFit(d, BPPARAM = BiocParallel::MulticoreParam(workers = 2))


### Fit null model proportions and test for sQTLs

# If possible, increase the number of workers
d <- dmTest(d, BPPARAM = BiocParallel::MulticoreParam(workers = 2))

plotTest(d)

dataSQTL_dmSQTLtest <- d

use_data(dataSQTL_dmSQTLtest, overwrite = TRUE)


### Plot feature proportions for top sQTL

res <- results(d)
res <- res[order(res$pvalue, decreasing = FALSE), ]

gene_id <- res$gene_id[1]
snp_id <- res$snp_id[1]


plotFit(d, gene_id, snp_id)











