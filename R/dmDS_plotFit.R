
# gene_id = gene_id; counts = x@counts; samples = x@samples; dispersion = slot(x, x@dispersion); proportions_full = x@fit_full; proportions_null = NULL; table = NULL; plot_type = "barplot"; order = FALSE; plot_full = TRUE; plot_null = FALSE; out_dir = "~/"


dmDS_plotFit <- function(gene_id, counts, samples, dispersion = numeric(), proportions_full = NULL, proportions_null = NULL, table = NULL, plot_type = c("barplot", "boxplot1", "boxplot2", "lineplot", "ribbonplot")[3], order = TRUE, plot_full = ifelse(is.null(proportions_full), FALSE, TRUE), plot_null = ifelse(is.null(proportions_null), FALSE, TRUE), plot_main = TRUE, out_dir = NULL){
  
  
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
    main <- NULL
    
    
    if(plot_main){
      
      mean_expression_gene <- mean(colSums(counts_gene), na.rm = TRUE)
      
      main <- paste0(gene, "\n Mean expression = ", round(mean_expression_gene))
      
      if(length(dispersion) > 0){
        
        if(length(dispersion) == 1)
          dispersion_gene <- dispersion
        else
          dispersion_gene <- dispersion[gene]
        
        main <- paste0(main, ", Dispersion = ", round(dispersion_gene, 2))
        
      }
      
      if(!is.null(table)){
        
        table_tmp <- table[table$gene_id == gene, ]
        
        main <- paste0(main, "\n LR = ", round(table_tmp["lr"], 2) , ", P-value = ", sprintf("%.02e", table_tmp["pvalue"]), ", FDR = ", sprintf("%.02e", table_tmp["adj_pvalue"]))    
        
      }
    }
    
    
    pi_full <- NULL
    pi_null <- NULL
    
    if(plot_full)
      pi_full <- proportions_full[[gene]]
    if(plot_null)
      pi_null <- proportions_null[[gene]]
    
    
    ggp <- dm_plotProportions(counts = counts_gene, group = group, pi_full = pi_full, pi_null = pi_null, main = main, plot_type = plot_type, order = order)
    
    
    
    if(!is.null(out_dir))
      pdf(paste0(out_dir, "dmfit_", gsub(pattern = "\\.", replacement = "_" , paste0(i, "_", gene)), ".pdf"), width = 12, height = 7)
    
    print(ggp)
    
    if(!is.null(out_dir))
      dev.off()
    
    
  }
  
  
  
}










