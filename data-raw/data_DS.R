# Prepare data for examples and vignette 

setwd("/home/gosia/R/multinomial_project/package_devel/DM/data-raw/simulations_sim5_drosophila_noDE_noNull/")

library(DM)

library(devtools)



########################################################
# load metadata
########################################################

# create metadata file
metadata <- data.frame(sample_id = paste0("sample_",1:6), group=c(rep("C1", 3), rep("C2", 3)), stringsAsFactors = FALSE)
metadata$group <- as.factor(metadata$group)

metadata



################################################################################
# htseq counts
################################################################################

htseq <- read.table("2_counts/dexseq_nomerge/htseq_counts.txt", header = TRUE, as.is = TRUE)
htseq <- htseq[!grepl(pattern = "_", htseq$group_id), ]


counts <- as.matrix(htseq[,-1])
group_id <- htseq[,1]
group_split <- limma::strsplit2(group_id, ":")
gene_id <- group_split[, 1]
feature_id <- group_split[, 2]
sample_id = metadata$sample_id
group = metadata$group




### sample only 100 random genes 
# set.seed(1)
# genes_subset <- unique(gene_id)[sample(length(unique(gene_id)), size = 100, replace = FALSE)]


d <- dmDSdata(counts = counts, gene_id = gene_id, feature_id = feature_id, sample_id = sample_id, group = group)

plotData(d)
dev.off()

d <- dmFilter(d, min_samps_gene_expr = 6)

plotData(d)
dev.off()


### sample only 100 random genes
set.seed(1)
genes_subset <- names(d)[sample(length(d), size = 100, replace = FALSE)]



# d <- dmDispersion(d, verbose = TRUE, BPPARAM = BiocParallel::MulticoreParam(workers = 20))

# plotDispersion(d)
# dev.off()

# genes_subset <- which(log10(d@genewise_dispersion) > 2 & log10(d@genewise_dispersion) < 6)
# length(genes_subset)



################################################################################
# data for DS examples
################################################################################



dataDS_counts <- htseq[gene_id %in% genes_subset, ]
dataDS_metadata <- metadata


use_data(dataDS_counts, dataDS_metadata, overwrite = TRUE)




### Start examples
counts <- as.matrix(dataDS_counts[,-1])
group_id <- dataDS_counts[,1]
group_split <- limma::strsplit2(group_id, ":")
gene_id <- group_split[, 1]
feature_id <- group_split[, 2]
sample_id = dataDS_metadata$sample_id
group = dataDS_metadata$group


d <- dmDSdata(counts = counts, gene_id = gene_id, feature_id = feature_id, sample_id = sample_id, group = group)

plotData(d)
dev.off()

### End examples


dataDS_dmDSdata <- d

use_data(dataDS_dmDSdata, overwrite = TRUE)




### Start examples

dd <- dataDS_dmDSdata
plot(dd)
plotData(dd)
dev.off()

d <- dmFilter(dd)
plot(dd)
plotData(d)
dev.off()

d <- dmFilter(dd, max_features = 10)
plotData(d)
dev.off()

### End examples




### Start examples

d <- dataDS_dmDSdata

d <- dmFilter(d)

# If possible, increase the number of workers
d <- dmDispersion(d, BPPARAM = BiocParallel::MulticoreParam(workers = 5))

### End examples



dataDS_dmDSdispersion <- d

use_data(dataDS_dmDSdispersion, overwrite = TRUE)



### Start examples

d <- dataDS_dmDSdispersion

plotDispersion(d)
plot(d)
dev.off()

### End examples



### Start examples
d <- dataDS_dmDSdispersion

# If possible, increase the number of workers
d <- dmFit(d, BPPARAM = BiocParallel::MulticoreParam(workers = 1))

gene_id <- names(d)[1]

plotFit(d, gene_id = gene_id)
plot(d, gene_id = gene_id, plot_type = "lineplot", plot_full = FALSE)
dev.off()

### End examples



### Start examples
d <- dataDS_dmDSdispersion

# If possible, increase the number of workers
d <- dmFit(d, BPPARAM = BiocParallel::MulticoreParam(workers = 1))

d <- dmLRT(d)

results <- results(d)

plotLRT(d)
plot(d)
dev.off()

gene_id <- results$gene_id[1]
plotFit(d, gene_id = gene_id)
dev.off()


### End examples





