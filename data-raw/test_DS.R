# Test DM on entire data set

setwd("/home/gosia/R/multinomial_project/package_devel/DM/data-raw/simulations_sim5_drosophila_noDE_noNull/")

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

htseq <- read.table("2_counts/dexseq_nomerge/htseq_counts.txt", header = TRUE, as.is = TRUE)
htseq <- htseq[!grepl(pattern = "_", htseq$group_id), ]


counts <- htseq[,-1]
group_id <- htseq[,1]
group_split <- limma::strsplit2(group_id, ":")
gene_id <- group_split[, 1]
feature_id <- group_split[, 2]
sample_id = metadata$sample_id
group = metadata$group



data <- dmDSdata(counts = counts, gene_id = gene_id, feature_id = feature_id, sample_id = sample_id, group = group)

plotData(data)
dev.off()



system.time(data <- dmFilter(data))

plotData(data)
dev.off()



data <- dmDispersion(data, verbose = TRUE, BPPARAM = BiocParallel::MulticoreParam(workers = 20))

plotDispersion(data)
dev.off()



data <- dmFit(data, BPPARAM = BiocParallel::MulticoreParam(workers = 15))


plotFit(data, gene_id = "FBgn0001316", plot_type = "barplot")
dev.off()


# save(data, file = "data.Rdata")
load("data.Rdata")


data <- dmLRT(data, BPPARAM = BiocParallel::MulticoreParam(workers = 10))

plotLRT(data, out_dir = "./")

results <- results(data)




plotFit(data, gene_id = results$gene_id[1:5], out_dir = "./")






data <- dmLRT(data, compared_groups = list(a = 1:2, b = 1:2), BPPARAM = BiocParallel::MulticoreParam(workers = 10))

plotLRT(data, out_dir = "./")

results <- results(data)



plotFit(data, gene_id = results$gene_id[1:5], out_dir = "./")





data <- dmLRT(data, compared_groups = list(a = 1:2), BPPARAM = BiocParallel::MulticoreParam(workers = 10))

plotLRT(data, out_dir = "./k")

results <- results(data)


plotFit(data, gene_id = results$gene_id[1:5], out_dir = "./k")








