# Prepare data for examples and vignette 

setwd("/home/gosia/R/multinomial_project/package_devel/DRIMSeq/")

library(DRIMSeq)

library(devtools)


## Input files: transcript expression, gene location and genotype information
data_dir <- "data-raw/geuvadis/data/"


########################################################
# process data downloaded from the GEUVADIS website 
########################################################

library(GenomicRanges)
library(rtracklayer)
library(limma)

### annotation

gtf0 <- import(paste0(data_dir, "geuvadis_annotation/gencode.v12.annotation.gtf"))

# keep protein coding genes
keep_index <- mcols(gtf0)$gene_type == "protein_coding" & mcols(gtf0)$type == "gene" 
gtf <- gtf0[keep_index]
# remove 'chr'
seqlevels(gtf) <- gsub(pattern = "chr", replacement = "", x = seqlevels(gtf))

genes_bed <- data.frame(chr = seqnames(gtf), start =  start(gtf), end = end(gtf), geneId = mcols(gtf)$gene_id)

write.table(genes_bed, paste0(data_dir, "annotation/genes.bed"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


### samples

samples <- read.table(paste0(data_dir, "geuvadis_analysis_results/E-GEUV-1.sdrf.txt"), header = T, sep="\t", as.is=TRUE)

table(samples$Assay.Name)

samples <- samples[c("Assay.Name", "Characteristics.population.")]
samples <- unique(samples)

dim(samples)
table(samples$Characteristics.population.)

colnames(samples) <- c("sample_id", "population")
samples$sample_id_short <- strsplit2(samples$sample_id, "\\.")[,1]



### expression in counts

expr_all <- read.table(paste0(data_dir, "geuvadis_analysis_results/GD660.TrQuantCount.txt"), header = T, sep="\t", as.is = TRUE)


table(colnames(expr_all) %in% samples$sample_id)


expr_all <- expr_all[, c("TargetID", "Gene_Symbol", samples$sample_id)]
colnames(expr_all) <- c("trId", "geneId", samples$sample_id_short)


for(i in unique(samples$population)){
  
  expr <- expr_all[, c("trId", "geneId", samples$sample_id_short[samples$population == i])]

  write.table(expr, paste0(data_dir, "expression/TrQuantCount_",i,".tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
}


### genotypes
library(Rsamtools)
library(VariantAnnotation)
library(tools)


files <- list.files(path = paste0(data_dir, "geuvadis_genotypes"), pattern = "genotypes.vcf.gz", full.names = TRUE, include.dirs = FALSE)

## bigzip and index the vcf files
for(i in 1:length(files)){
  # i = 1
  
  zipped <- bgzip(files[i])
  idx <- indexTabix(zipped, format = "vcf")
  
}

## extended gene ranges
window <- 5000
gene_ranges <- resize(gtf, GenomicRanges::width(gtf) + 2 * window, fix = "center")

population <- unique(samples$population)
chr <- gsub("chr", "", strsplit2(files, split = "\\.")[, 2])


for(j in 1:length(population)){
  
  for(i in 1:length(files)){
    
    cat(population[j], chr[i], fill = TRUE)
    
    zipped <- paste0(file_path_sans_ext(files[i]), ".bgz")
    idx <- paste0(file_path_sans_ext(files[i]), ".bgz.tbi")
    tab <- TabixFile(zipped, idx)
    
    ## Explore the file header with scanVcfHeader
    hdr <- scanVcfHeader(tab)
    print(all(samples$sample_id_short %in% samples(hdr)))
    
    ## Read VCF file 
    gene_ranges_tmp <- gene_ranges[seqnames(gene_ranges) == chr[i]]
    
    param <- ScanVcfParam(which = gene_ranges_tmp, samples = samples$sample_id_short[samples$population == population[j]])
    vcf <- readVcf(tab, "hg19", param)
    
    
    ## Keep only the bi-allelic SNPs
    
    # width of ref seq
    rw <- width(ref(vcf))
    # width of first alt seq
    aw <- unlist(lapply(alt(vcf), function(x) {width(x[1])}))
    # number of alternate genotypes
    nalt <- elementLengths(alt(vcf))
    # select only bi-allelic SNPs (monomorphic OK, so aw can be 0 or 1)
    snp <- rw == 1 & aw <= 1 & nalt == 1
    # subset vcf  
    vcfbi <- vcf[snp,]
    
    rowdata <- rowData(vcfbi)
    
    ## Convert genotype into number of alleles different from reference
    geno <- geno(vcfbi)$GT
    geno01 <- geno
    geno01[,] <- -1
    geno01[geno %in% c("0/0", "0|0")] <- 0 # REF/REF
    geno01[geno %in% c("0/1", "0|1", "1/0", "1|0")] <- 1 # REF/ALT
    geno01[geno %in% c("1/1", "1|1")] <- 2 # ALT/ALT
    # geno01 should be integer, not character
    mode(geno01) <- "integer"
    
    genotypes <- unique(data.frame(chr = seqnames(rowdata), start = start(rowdata), end = end(rowdata), snpId = rownames(geno01), geno01, stringsAsFactors = FALSE))
    
    ### sorting
    genotypes <- genotypes[order(genotypes[,2]), ]
    
    write.table(genotypes, file = paste0(data_dir, "genotypes/genotypes_", j, "_chr", chr[i], ".tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    
    
  }
}


########################################################
# Choose a subset of genes 
########################################################


data_dir <- "data-raw/geuvadis/data/"

### Create dmSQTLdata object


# gene_ranges with names!
gene_ranges <- rtracklayer::import(paste0(data_dir, "annotation/genes.bed"))
names(gene_ranges) <- S4Vectors::mcols(gene_ranges)$name

counts <- read.table(paste0(data_dir, "expression/TrQuantCount_CEU.tsv"), header = TRUE, sep = "\t", as.is = TRUE)

genotypes <- read.table(paste0(data_dir, "genotypes/genotypes_CEU_chr19.tsv"), header = TRUE, sep = "\t", as.is = TRUE)

# snp_ranges with names!
snp_ranges <- GenomicRanges::GRanges(S4Vectors::Rle(genotypes$chr), IRanges::IRanges(genotypes$start, genotypes$end))
names(snp_ranges) <- genotypes$snpId 

## Check if samples in count and genotypes are in the same order
all(colnames(counts[, -(1:2)]) == colnames(genotypes[, -(1:4)]))
sample_id <- colnames(counts[, -(1:2)])


d <- d_org <- dmSQTLdataFromRanges(counts = counts[, -(1:2)], gene_id = counts$geneId, feature_id = counts$trId, gene_ranges = gene_ranges, genotypes = genotypes[, -(1:4)], snp_id = genotypes$snpId, snp_ranges = snp_ranges, sample_id = sample_id, window = 5e3, BPPARAM = BiocParallel::MulticoreParam(workers = 2))

plotData(d, out_dir = "./")



### Filtering
d <- dmFilter(d, min_samps_gene_expr = 70, min_samps_feature_prop = 5, minor_allele_freq = 5, BPPARAM = BiocParallel::MulticoreParam(workers = 2))

plotData(d, out_dir = "./")


# set.seed(123)
# genes_subset <- names(d)[sample(length(d), size = 50, replace = FALSE)]

oo <- order(width(d@genotypes), decreasing = FALSE)
genes_subset <- names(d)[oo][1:50]



d_sub1 <- d_org[genes_subset, ]
snps_subset1 <- unique(d_sub1@blocks@unlistData[, "snp_id"])


d_sub2 <- d[genes_subset, ]
snps_subset2 <- unique(d_sub2@blocks@unlistData[, "snp_id"])


snps_extra <- setdiff(snps_subset1, snps_subset2)
snps_extra <- snps_extra[sample(length(snps_extra), size = 500, replace = FALSE)]


write.table(genes_subset, paste0("inst/extdata/gene_id_subset.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(c(snps_subset2, snps_extra), paste0("inst/extdata/snp_id_subset.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


########################################################
# Keep data only for 50 genes from chr19 
########################################################

# data_ex_dir <- system.file("extdata", package = "DRIMSeq")

# genes_subset = readLines(file.path(data_ex_dir, "/gene_id_subset.txt"))
# snps_subset = readLines(file.path(data_ex_dir, "/snp_id_subset.txt"))


genes_subset = readLines(file.path("inst/extdata/gene_id_subset.txt"))
snps_subset = readLines(file.path("inst/extdata/snp_id_subset.txt"))


### Subset gene ranges
gene_bed <- read.table(paste0(data_dir, "annotation/genes.bed"), header = FALSE, as.is = TRUE)
gene_bed <- gene_bed[gene_bed[, 4] %in% genes_subset, ]

write.table(gene_bed, paste0("inst/extdata/genes_subset.bed"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


### Subset counts
counts <- read.table(paste0(data_dir, "expression/TrQuantCount_CEU.tsv"), header = TRUE, as.is = TRUE)
counts <- counts[counts$geneId %in% genes_subset, ]

write.table(counts, paste0("inst/extdata/TrQuantCount_CEU_subset.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


### Subset genotypes
genotypes <- read.table(paste0(data_dir, "genotypes/genotypes_CEU_chr19.tsv"), header = TRUE, sep = "\t", as.is = TRUE)
genotypes <- genotypes[genotypes$snpId %in% snps_subset, ]

write.table(genotypes, paste0("inst/extdata/genotypes_CEU_subset.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


########################################################
# Examples
########################################################


#############################
### Create dmSQTLdata object
#############################

library(GenomicRanges)
library(rtracklayer)

data_dir  <- system.file("extdata", package = "DRIMSeq")


# gene_ranges with names!
gene_ranges <- import(paste0(data_dir, "/genes_subset.bed"))
names(gene_ranges) <- mcols(gene_ranges)$name

counts <- read.table(paste0(data_dir, "/TrQuantCount_CEU_subset.tsv"), header = TRUE, sep = "\t", as.is = TRUE)

genotypes <- read.table(paste0(data_dir, "/genotypes_CEU_subset.tsv"), header = TRUE, sep = "\t", as.is = TRUE)

# snp_ranges with names!
snp_ranges <- GRanges(Rle(genotypes$chr), IRanges(genotypes$start, genotypes$end))
names(snp_ranges) <- genotypes$snpId 

## Check if samples in count and genotypes are in the same order
all(colnames(counts[, -(1:2)]) == colnames(genotypes[, -(1:4)]))
sample_id <- colnames(counts[, -(1:2)])


d <- dmSQTLdataFromRanges(counts = counts[, -(1:2)], gene_id = counts$geneId, feature_id = counts$trId, gene_ranges = gene_ranges, genotypes = genotypes[, -(1:4)], snp_id = genotypes$snpId, snp_ranges = snp_ranges, sample_id = sample_id, window = 5e3, BPPARAM = BiocParallel::MulticoreParam(workers = 2))

plotData(d)


data_dmSQTLdata <- d
use_data(data_dmSQTLdata, overwrite = TRUE)

#############################
### sQTL analysis
#############################
# If possible, increase the number of workers in BPPARAM

d <- data_dmSQTLdata

head(names(d))
length(d)
d[1:10, ]
d[1:10, 1:10]

### Filtering
d <- dmFilter(d, min_samps_gene_expr = 70, min_samps_feature_prop = 5, minor_allele_freq = 5, BPPARAM = BiocParallel::MulticoreParam(workers = 2))
plotData(d)


### Calculate dispersion
d <- dmDispersion(d, BPPARAM = BiocParallel::MulticoreParam(workers = 2))
plotDispersion(d)


### Fit full model proportions
d <- dmFit(d, BPPARAM = BiocParallel::MulticoreParam(workers = 2))


### Fit null model proportions and test for sQTLs
d <- dmTest(d, BPPARAM = BiocParallel::MulticoreParam(workers = 2))
plotTest(d)

head(results(d))

### Plot feature proportions for top sQTL
res <- results(d)
res <- res[order(res$pvalue, decreasing = FALSE), ]

gene_id <- res$gene_id[1]
snp_id <- res$snp_id[1]

plotFit(d, gene_id, snp_id)
plotFit(d, gene_id, snp_id, plot_type = "boxplot2", order = FALSE)
plotFit(d, gene_id, snp_id, plot_type = "ribbonplot")










