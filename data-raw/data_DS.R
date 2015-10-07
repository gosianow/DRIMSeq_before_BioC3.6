# Prepare data for examples and vignette 

# setwd("/home/gosia/R/multinomial_project/package_devel/DM/data-raw/simulations_sim5_drosophila_noDE_noNull/")


library(DM)

library(devtools)

data_path <- "data-raw/simulations_sim5_drosophila_noDE_noNull/"

########################################################
# load metadata
########################################################

# create metadata file
metadata <- data.frame(sample_id = paste0("sample_",1:6), group=c(rep("C1", 3), rep("C2", 3)), stringsAsFactors = FALSE)
metadata$group <- as.factor(metadata$group)

metadata

########################################################
# Simulation details
########################################################


simulation_details <- read.table(paste0(data_path, "3_truth/simulation_details.txt"), header = TRUE, as.is = TRUE)

status <- unique(simulation_details[, c("gene_id", "gene_ds_status")])

colnames(status) <- c("gene_id", "status")

table(status$status)

################################################################################
# htseq counts
################################################################################

htseq <- read.table(paste0(data_path, "2_counts/dexseq_nomerge/htseq_counts.txt"), header = TRUE, as.is = TRUE)
htseq <- htseq[!grepl(pattern = "_", htseq$group_id), ]


counts <- as.matrix(htseq[,-1])
group_id <- htseq[,1]
group_split <- limma::strsplit2(group_id, ":")
gene_id <- group_split[, 1]
feature_id <- group_split[, 2]
sample_id = metadata$sample_id
group = metadata$group

h(Encoding(gene_id))

### sample only 100 random genes 

d <- dmDSdata(counts = counts, gene_id = gene_id, feature_id = feature_id, sample_id = sample_id, group = group)

plotData(d)
dev.off()

d <- dmFilter(d, min_samps_gene_expr = 6)

plotData(d)
dev.off()


### sample only 100 random genes: 80 with status 0 and 20 with status 1
set.seed(1)

genes_subset0 <- status[status$status == 0, "gene_id"][sample(sum(status$status == 0), size = 80, replace = FALSE)]

set.seed(1)

genes_subset1 <- status[status$status == 1, "gene_id"][sample(sum(status$status == 1), size = 20, replace = FALSE)]


genes_subset <- c(genes_subset0, genes_subset1)


# d <- dmDispersion(d, verbose = TRUE, BPPARAM = BiocParallel::MulticoreParam(workers = 3))

# plotDispersion(d)
# dev.off()

# genes_subset <- which(log10(d@genewise_dispersion) > 2 & log10(d@genewise_dispersion) < 6)
# length(genes_subset)



################################################################################
# data for DS examples
################################################################################



dataDS_counts <- htseq[gene_id %in% genes_subset, ]
rownames(dataDS_counts) <- NULL
dataDS_metadata <- metadata
dataDS_status <- status[status$gene_id %in% genes_subset, ]
rownames(dataDS_status) <- NULL

use_data(dataDS_counts, dataDS_metadata, dataDS_status, overwrite = TRUE)




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
dev.off()

### End examples



### Start examples
d <- dataDS_dmDSdispersion

# If possible, increase the number of workers
d <- dmFit(d, BPPARAM = BiocParallel::MulticoreParam(workers = 2))

gene_id <- names(d)[1]

plotFit(d, gene_id = gene_id)
dev.off()

### End examples


dataDS_dmDSfit <- d

use_data(dataDS_dmDSfit, overwrite = TRUE)




### Start examples
d <- dataDS_dmDSfit

d <- dmTest(d, BPPARAM = BiocParallel::MulticoreParam(workers = 2))

res <- results(d)
res <- res[order(res$pvalue, decreasing = FALSE), ]

plotTest(d)
dev.off()

gene_id <- res$gene_id[1]
plotFit(d, gene_id = gene_id)
dev.off()


### End examples



dataDS_dmDStest <- d

use_data(dataDS_dmDStest, overwrite = TRUE)



status <- dataDS_status

d <- dataDS_dmDStest
results <- list()
results[[1]] <- results(d)

metadata <- data.frame(method = "DM")

data_ROCx <- calculateROCx(results, status)


plot_var = "method"
facet_var = NULL
plot_colors = NULL
xylim_one = TRUE

plotROCx(data_ROCx, metadata, plot_var = "method", facet_var = NULL, plot_colors = NULL, xylim_one = TRUE)



data_TPRFDR <- calculateTPRFDR(results, status, thresholds = c(0.01, 0.05, 0.1))

plotTPRFDR(data_TPRFDR, metadata, plot_var = "method", facet_var = NULL, plot_colors = NULL, xylim_one = TRUE)


data_venn <- calculateVenn(results, status = status, threshold = 0.05)


plotVenn(data_venn, plot_results = 1, metadata, plot_var = "method", plot_colors = NULL, plot_status = TRUE)








