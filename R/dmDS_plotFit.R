
# plot_type = c("barplot", "boxplot1", "boxplot2", "lineplot", "ribbonplot")[3]; order = TRUE; plot_full = TRUE; plot_null = TRUE; out_dir = "./"

dmDS_plotFit <- function(gene_id, counts, samples, dispersion = numeric(), proportions_full = NULL, proportions_null = NULL, table = NULL, plot_type = c("barplot", "boxplot1", "boxplot2", "lineplot", "ribbonplot")[3], order = TRUE, plot_full = ifelse(is.null(proportions_full), FALSE, TRUE), plot_null = ifelse(is.null(proportions_null), FALSE, TRUE), out_dir = NULL){


  for(i in 1:length(gene_id)){
    # i = 1
    cat(paste0("Plot gene ", i, ": ", gene_id[i], "\n"))
    
    gene <- gene_id[i]
    counts_gene <- counts[[gene]]
    
    if(nrow(counts_gene) <= 1){
      cat(paste0("!Gene has to have at least 2 features! \n"))
      next
    }
    
    group <- samples$group
    sample_id <- samples$sample_id
    mean_expression_gene <- mean(colSums(counts_gene), na.rm = TRUE)
    
    main <- paste0(gene, "\n Mean expression = ", round(mean_expression_gene))
    
    if(length(dispersion) > 0){
      
      if(length(dispersion) == 1)
      dispersion_gene <- dispersion
      else
      dispersion_gene <- dispersion[gene]
      
      main <- paste0(main, " / Dispersion = ", round(dispersion_gene, 2))
      
    }
    
    if(!is.null(table)){
      
      main <- paste0(main, "\n LR = ", round(table[gene, "lr"], 2) , " / P-value = ", sprintf("%.02e", table[gene, "pvalue"]), " / FDR = ", sprintf("%.02e", table[gene, "adj_pvalue"]))    
      
    }
    
    pi_full <- NULL
    pi_null <- NULL
    
    if(plot_full)
    pi_full <- proportions_full[[gene]]
    if(plot_null)
    pi_null <- proportions_null[[gene]]


    ggp <- dm_plotProportion(counts = counts_gene, group = group, sample_id = sample_id, pi_full = pi_full, pi_null = pi_null, main = main, plot_type = plot_type, order = order)



    if(!is.null(out_dir))
    pdf(paste0(out_dir, "proportions_", gsub(pattern = "\\.", replacement = "_" , gene), ".pdf"), width = 12, height = 7)
    
    print(ggp)
    
    if(!is.null(out_dir))
    dev.off()


  }

}










