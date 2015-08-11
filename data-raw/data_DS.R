# Prepare data for examples and vignette 


# library(devtools); library(GenomicRanges); library(BiocParallel); library(edgeR)


# Rfiles <- list.files("R/", full.names=TRUE); for(i in Rfiles) source(i)


########################################################
# load metadata
########################################################

# create metadata file
metadata <- data.frame(sample_id = paste0("sample_",1:6), group=c(rep("C1", 3), rep("C2", 3)), stringsAsFactors = FALSE)
metadata$group <- as.factor(metadata$group)

metadata



##############################################################################################################
# htseq counts
##############################################################################################################

htseq <- read.table("data-raw/simulations_sim5_drosophila_noDE_noNull/2_counts/dexseq_nomerge/htseq_counts.txt", header = TRUE, as.is = TRUE)
htseq <- htseq[!grepl(pattern = "_", htseq$group_id), ]


counts <- as.matrix(htseq[,-1])
group_id <- htseq[,1]
group_split <- limma::strsplit2(group_id, ":")
gene_id_counts <- group_split[, 1]
feature_id_counts <- group_split[, 2]
sample_id = metadata$sample_id
group = metadata$group


### sample only 100 genes 

length(unique(gene_id_counts))

set.seed(1)
genes_subset <- unique(gene_id_counts)[sample(length(unique(gene_id_counts)), size = 100, replace = FALSE)]


dataDS_counts <- htseq[gene_id_counts %in% genes_subset, ]
dataDS_metadata <- metadata

use_data(dataDS_counts, dataDS_metadata, overwrite = TRUE)




counts <- as.matrix(dataDS_counts[,-1])
group_id <- dataDS_counts[,1]
group_split <- limma::strsplit2(group_id, ":")
gene_id_counts <- group_split[, 1]
feature_id_counts <- group_split[, 2]
sample_id = dataDS_metadata$sample_id
group = dataDS_metadata$group


data_org <- DM::dmDSdata(counts = counts, gene_id_counts = gene_id_counts, feature_id_counts = feature_id_counts, sample_id = sample_id, group = group)

dmDSplotData(data_org)


dataDS_dmDSdata <- data_org

use_data(dataDS_dmDSdata, overwrite = TRUE)


######################################
# htseq counts: filtering
######################################


data <- dataDS_dmDSdata
dmDSplotData(data)


data <- dmDSfilter(data)
dmDSplotData(data)


data <- dmDSdispersion(data)


dmDSplotDispersion(data)


data_fit <- dmDSfit(data, dispersion = "tagwise_dispersion", prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = MulticoreParam(workers = 1))



dmDSplotFit(data_fit, gene_id =   , plot_type = "barplot", order = TRUE, plot_full = FALSE, plot_nunll = FALSE)



table <- dmDStest(data_fit)


dmDSplotTest(table)


