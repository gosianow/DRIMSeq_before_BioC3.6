
#' Venn diagram
#'  
#' @inheritParams calculateROCx
#' @param threshold Threshold for the FDR.
#' @return 
#' calculateVenn returns a list of vectors with names of the genes that are significant \code{adj_pvalue < threshold}. If status is not \code{NULL}, the last element of the output is a list of genes with status equal 1.
#' 
#' @seealso \code{\link{plotROCx}}, \code{\link{plotTPRFDR}}
#'  @author Malgorzata Nowicka
#' @export
calculateVenn <- function(results, status = NULL, threshold = 0.05){
  
  venn <- lapply(results, function(r){
    
    if(!all(c("gene_id", "adj_pvalue") %in% colnames(r))){
      
      message("Some columns 'gene_id', 'adj_pvalue' are missing in one of the results")
      character()
      
    }else{
      
      as.character(r[which(r$adj_pvalue < threshold), "gene_id"])
      
    }

    })
  
  names(venn) <- NULL
  
  if(!is.null(status)){
    stopifnot(all(c("gene_id", "status") %in% colnames(status)))
    venn[["status"]] <- as.character(status[which(status$status == 1), "gene_id"])
  }

  
  return(venn)
}



################################################################################

#' @param data_venn List structured like the output of \code{calculateVenn}.
#' @param plot_var Name of one of the variables in \code{metadata}. Levels of this variable will be used as category names in venn diagram.
#' @param plot_results Numeric vector of maximum length 4 (if \code{plot_status = TRUE}) or 5 (if \code{plot_status = FALSE}) indicating which results should be plotted. 
#' @param plot_status Logical. Whether to plot a circle that corresponds to truly significant genes.
#' @examples 
#' 
#' status <- dataDS_status
#' 
#' d <- dataDS_dmDStest
#' results <- list()
#' results[[1]] <- results(d)
#' 
#' metadata <- data.frame(method = "DM")
#' 
#' data_venn <- calculateVenn(results, status = status, threshold = 0.05)
#' 
#' plotVenn(data_venn, plot_results = 1, metadata, plot_var = "method", 
#' plot_colors = NULL, plot_status = TRUE)
#' 
#' @rdname calculateVenn
#' @export
plotVenn <- function(data_venn, plot_results, metadata, plot_var, plot_colors = NULL, plot_status = TRUE){
  
  stopifnot(is.logical(plot_status))
  if("status" %in% names(data_venn))
    stopifnot(length(data_venn) == nrow(metadata) + 1)
  else
    stopifnot(length(data_venn) == nrow(metadata))
  stopifnot(plot_var %in% colnames(metadata))
  
  if(!is.null(plot_colors))
    stopifnot(nlevels(metadata[, plot_var]) == length(plot_colors))
  
  if(plot_status && "status" %in% names(data_venn)){
    stopifnot(length(plot_results) <= 4)
    
    data_venn_tmp <- data_venn[plot_results]
    names(data_venn_tmp) <- metadata[plot_results, plot_var]
    data_venn_tmp[["status"]] <- data_venn[["status"]]
    
    if(!is.null(plot_colors)){
      colors <- metadata[, plot_var]
      levels(colors) <- plot_colors
      colors <- c(colors[plot_results], "grey")
    }else{
      colors <- c(colorb(length(plot_results)), "grey")
    }

  }else{
    stopifnot(length(plot_results) <= 5)
    
    data_venn_tmp <- data_venn[plot_results]
    names(data_venn_tmp) <- metadata[plot_results, plot_var]
    
    if(!is.null(plot_colors)){
      colors <- metadata[, plot_var]
      levels(colors) <- plot_colors
      colors <- c(colors[plot_results])
    }else{
      colors <- c(colorb(length(plot_results)))
    }
    
  }

  cex = 1.5
  cat.cex=1.5
  lwd=2
  lty=1
  alpha=0.5
  margin=0.1
  
  if(length(colors) == 4)
    colors <- colors[c(1, 3, 4 ,2)]
  
  venn <- VennDiagram::venn.diagram(data_venn_tmp, filename = NULL, fill = colors, col = colors, cex=cex, cat.cex = cat.cex, lwd=lwd, lty=lty, alpha=alpha, margin = margin, scaled = FALSE)
  
  
  return(grid::grid.draw(venn))
  
  
}

