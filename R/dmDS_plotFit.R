
# plot_type = c("barplot", "boxplot1", "boxplot2", "lineplot")[2]; order = TRUE; fit = NULL; table = NULL; plot_full = ifelse(is.null(fit), FALSE, TRUE); plot_null = ifelse(is.null(fit), FALSE, TRUE); out_dir = "./"

dmDS_plotFit <- function(data, gene_id, plot_type = c("barplot", "boxplot1", "boxplot2", "lineplot")[2], order = TRUE, fit = NULL, table = NULL, plot_full = ifelse(is.null(fit), FALSE, TRUE), plot_null = ifelse(is.null(fit), FALSE, TRUE), out_dir = "./"){


  for(i in 1:length(gene_id)){
    # i = 1
    cat(paste0("Plot gene ", i, ": ", gene_id[i], "\n"))
    
    gene <- gene_id[i]
    counts <- data@counts[[gene]]
    group <- data@samples$group
    sample_id <- data@samples$sample_id
    mean_expression_gene <- mean(colSums(counts), na.rm = TRUE)
    
fit_full <- fit_null <- NULL

    if(is.null(fit)){      
        dispersion_gene <- NA
        }else{
            dispersion_gene <- fit$fit_full[[gene]]$gamma0
            if(plot_full)
            fit_full <- fit$fit_full[[gene]]
            if(plot_null)
            fit_null <- fit$fit_null[[gene]]          
        }


        if(is.null(table)){
          main <- paste0(gene, "\n Mean expression = ", round(mean_expression_gene), " / Dispersion = ", round(dispersion_gene, 2))
          }else{
            main <- paste0(gene, "\n Mean expression = ", round(mean_expression_gene), " / Dispersion = ", round(dispersion_gene, 2), "\n LR = ", round(table[gene, "lr"], 2) , " / P-value = ", sprintf("%.02e", table[gene, "pvalue"]), " / FDR = ", sprintf("%.02e", table[gene, "adj_pvalue"]))

        }


        ggp <- dm_plotProportion(counts, group, sample_id, fit_full = fit_full, fit_null = fit_null, main = main, plot_type = plot_type, order = order)



        pdf(paste0(out_dir, "proportions_", gsub(pattern = "\\.", replacement = "_" , gene_id), ".pdf"), width = 12, height = 7)
        print(ggp)
        dev.off()


    }





}










