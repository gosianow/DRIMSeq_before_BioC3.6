
# plot_type = c("barplot", "boxplot1", "boxplot2", "lineplot", "ribbonplot")[3]; order = TRUE; plot_full = TRUE; plot_null = TRUE; out_dir = "./"

dmSQTL_plotFit <- function(data, gene_id, snp_id, fit = NULL, table = NULL, plot_type = c("barplot", "boxplot1", "boxplot2", "lineplot", "ribbonplot")[3], order = TRUE, plot_full = ifelse(is.null(fit), FALSE, TRUE), plot_null = ifelse(is.null(fit), FALSE, TRUE), out_dir = "./"){


  for(i in 1:length(gene_id)){
    # i = 1
    cat(paste0("Plot pair ", i, ": ", gene_id[i], ":", snp_id[i], "\n"))
    
    gene <- gene_id[i]
    snp <- snp_id[i]
    counts <- data@counts[[gene]]
    group <- data@genotypes[[gene]][snp, ]

    NAs <- !(is.na(counts[1,]) | is.na(group))
    counts <- counts[, NAs, drop = FALSE]
    group <- factor(group[NAs])
    sample_id <- data@samples$sample_id[NAs]

    mean_expression_gene <- mean(colSums(counts), na.rm = TRUE)
    
    fit_full <- fit_null <- NULL

    if(is.null(fit)){      
        dispersion_gene <- NA
        }else{
            dispersion_gene <- fit$fit_full[[gene]][[snp]]$gamma0
            if(plot_full)
            fit_full <- fit$fit_full[[gene]][[snp]]
            if(plot_null)
            fit_null <- fit$fit_null[[gene]][[snp]]        
        }



        if(is.null(table)){
          main <- paste0(gene,":", snp, "\n Mean expression = ", round(mean_expression_gene), " / Dispersion = ", round(dispersion_gene, 2))
          }else{
            index_table <- paste0(gene, ":", snp)
            main <- paste0(gene, ":", snp,"\n Mean expression = ", round(mean_expression_gene), " / Dispersion = ", round(dispersion_gene, 2), "\n LR = ", round(table[index_table, "lr"], 2) , " / P-value = ", sprintf("%.02e", table[index_table, "pvalue"]), " / FDR = ", sprintf("%.02e", table[index_table, "adj_pvalue"]))

        }


        ggp <- dm_plotProportion(counts, group, sample_id, fit_full = fit_full, fit_null = fit_null, main = main, plot_type = plot_type, order = order)


        pdf(paste0(out_dir, "proportions_", gsub(pattern = "\\.", replacement = "_" , paste0(gene, "_", snp)), ".pdf"), width = 12, height = 7)
        print(ggp)
        dev.off()


    }





}










