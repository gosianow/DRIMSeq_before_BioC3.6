
#' @export
calculate_venn <- function(results, status, threshold = 0.05){
  
  
  venn <- lapply(results, function(r){
    
    r[which(r$adj_pvalue < threshold), "gene_id"]
    
    })
  
  venn[["status"]] <- status[which(status$status == 1), "gene_id"]
  
  return(venn)
}



################################################################################

#' @export
plot_venn <- function(data_venn, plot_results, split_levels, plot_levels, plot_colors = NULL, plot_status = TRUE){
  
  
  if(plot_status && "status" %in% names(data_venn)){
    stopifnot(length(plot_results) <= 4)
    
    data_venn_tmp <- data_venn[plot_results]
    names(data_venn_tmp) <- split_levels[plot_results, plot_levels]
    data_venn_tmp[["status"]] <- data_venn[["status"]]
    
    if(!is.null(plot_colors))
    colors <- c(plot_colors, "grey")
    else
    colors <- c(colorb(length(plot_results)), "grey")
    
  }else{
    stopifnot(length(plot_results) <= 5)
    
    data_venn_tmp <- data_venn[plot_results]
    names(data_venn_tmp) <- split_levels[plot_results, plot_levels]
    
    if(!is.null(plot_colors))
    colors <- plot_colors
    else
    colors <- colorb(length(plot_results))
  }

  
  
  cex = 1.5
  cat.cex=1.5
  lwd=2
  lty=1
  alpha=0.5
  margin=0.1
  
  if(length(colors) == 4)
    colors <- colors[c(1, 3, 4 ,2)]
  
  venn <- VennDiagram::venn.diagram(data_venn_tmp, filename = NULL, fill = colors, col = colors, cex=cex, cat.cex=cat.cex, lwd=lwd, lty=lty, alpha=alpha, margin = margin, scaled = FALSE)
  
  
  return(grid::grid.draw(venn))
  
  
}

