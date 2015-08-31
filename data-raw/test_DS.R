# Test DM on entire data set

setwd("/home/gosia/R/multinomial_project/package_devel/DM")

library(DM)

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



data <- dmDSdata(counts = counts, gene_id_counts = gene_id_counts, feature_id_counts = feature_id_counts, sample_id = sample_id, group = group)

# dmDSplotData(data)


system.time(data_filt <- dmDSfilter(data))

# system.time(data_filt <- myFilter(data@counts,data@samples))



dmDSplotData(data)


data <- dmDSdispersion(data, BPPARAM = BiocParallel::MulticoreParam(workers = 3))

dmDSplotDispersion(data)


data <- dmDSfit(data, BPPARAM = BiocParallel::MulticoreParam(workers = 3))

dmDSplotFit(data, gene_id = "FBgn0001316", plot_type = "barplot")


data <- dmDStest(data)

dmDSplotTest(data)







