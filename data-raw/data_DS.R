# Prepare data for examples and vignette 



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



######################################
# htseq counts: filtering
######################################


min_samps_feature_prop_list <- 3
min_feature_prop_list <- 0.01 # in cpm
min_samps_gene_expr_list <- 6
min_gene_expr_list <- 1
max_features_list <- Inf


i = 1

min_samps_gene_expr = min_samps_gene_expr_list[i]
min_gene_expr = min_gene_expr_list[i]
min_samps_feature_prop = min_samps_feature_prop_list[i]
min_feature_prop = min_feature_prop_list[i]
max_features = max_features_list[i]


data <- dmDSfilter(data_org, min_samps_gene_expr = min_samps_gene_expr, min_gene_expr = min_gene_expr, min_samps_feature_prop = min_samps_feature_prop, min_feature_prop = min_feature_prop, max_features = max_features)


dmDSplotData(data)

dataDS_dmDSdata <- data

use_data(dataDS_dmDSdata, overwrite = TRUE)






dispersion <- dmDSdispersion(data, mean_expression = TRUE, common_dispersion = TRUE, tagwise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+5), disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "none", disp_prior_df = 10, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = MulticoreParam(workers = 2))


dmDSplotDispersion(dispersion)


data_fit <- dmDSfit(data, dispersion = "tagwise_dispersion", prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = MulticoreParam(workers = 1))



dmDSplotFit(data_fit, gene_id =   , plot_type = "barplot", order = TRUE, plot_full = FALSE, plot_nunll = FALSE)



table <- dmDStest(data_fit)


dmDSplotTest(table)


