#' ROC plot
#' 
#' Creates ROC curves for a list of result data frames.
#' 
#' @param results List of data frames with results from different methods. Each data frame should contain \code{gene_id}, \code{pvalue} and \code{adj_pvalue} columns. Otherwise NAs are returned.
#' @param status Data frame with truth. It must have \code{gene_id} and \code{status} columns, where \code{status} can take values of 0 and 1 (or FALSE, TRUE).
#' @return 
#' \code{calculateROCx} returns a list where each element contains a data frame ROC with TPR and FPR and a data frame X with TPR and FPR calculated when the estimated FDR is lower than 0.05. 
#' 
#' \code{plotROCx} returns a plot with ROC curves. Additionally, on each curve, the TPR for the estimated FDR lower than 0.05 is marked with X.
#' 
#'  @author Malgorzata Nowicka
#'  @name plotROCx
#'  @export
calculateROCx <- function(results, status){
  
  stopifnot(all(c("gene_id", "status") %in% colnames(status)))
  
  status <- status[complete.cases(status[, c("gene_id", "status")]), , drop = FALSE]
  
  P <- sum(status$status == 1)
  N <- sum(status$status == 0)
  
  
  ROCXlist <- lapply(1:length(results), function(m){
    # m = 2
    # print(m)
    
    if(!all(c("gene_id", "pvalue", "adj_pvalue") %in% colnames(results[[m]]))){
      
      message("Some columns 'gene_id', 'pvalue', 'adj_pvalue' are missing in results ", m, "\nReturn NAs")
      
      ROC <- data.frame(FPR = NA, TPR = NA)
      X <- data.frame(FPR = NA, TPR = NA) 
      
    }else{
      
      sc <- merge(status, results[[m]], by = "gene_id", all.x = TRUE)
      
      NAs <- is.na(sc[, "pvalue"]) | is.na(sc[, "adj_pvalue"])
      
      pvs <- sc$pvalue[!NAs]
      apvs <- sc$adj_pvalue[!NAs]
      sts <- sc$status[!NAs]
      
      ord <- order(pvs, decreasing = FALSE)
      sts <- sts[ord]
      pvs <- pvs[ord]
      apvs <- apvs[ord]
      
      TPRv <- cumsum(sts) / P
      FPRv <- cumsum(!sts) / N
      ROC <- data.frame(FPR = FPRv, TPR = TPRv)
      
      ### what is the TPR if the error was controlled 
      TPR <- sum(apvs[sts == 1] < 0.05) / P    
      X <- data.frame(FPR = approx(TPRv, FPRv, xout=TPR)$y, TPR = TPR)    
      
    }
    
    return(list(ROC = ROC, X = X))
    
  })
  
  
  return(ROCXlist)
  
}


################################################################################

#' @param data_ROCx List structured like the output of \code{calculateROCx}.
#' @param metadata Data frame where each row corresponds to one element of \code{data_ROCx} list. Columns describe the origin of this results, for example, which method was used to obtain them.
#' @param plot_var Name of one of the variables in \code{metadata} that indicates which results should be displayed in a single panel. 
#' @param facet_var Vector of one or two variables in \code{metadata} which indicates how the results should be splitted into different panels. If \code{NULL}, no splitting is done.
#' @param plot_colors Vector of colors corresponding to the levels of variables indicated by \code{plot_var}. If \code{NULL}, default colors are used. 
#' @param xylim_one Logical. Whether to set the x and y limits to (0,1).
#' 
#' @examples 
#' 
#' status <- dataDS_status
#' d <- dataDS_dmDStest
#' 
#' results <- list()
#' results[[1]] <- results(d)
#' metadata <- data.frame(method = "DM")
#' 
#' data_ROCx <- calculateROCx(results, status)
#' 
#' plotROCx(data_ROCx, metadata, plot_var = "method", facet_var = NULL, 
#'    plot_colors = NULL, xylim_one = TRUE)
#' 
#' @seealso \code{\link{dataDS_status}}, \code{\link{dataDS_dmDStest}}, \code{\link{plotTPRFDR}}, \code{\link{plotVenn}}
#' @rdname plotROCx
#' @export
plotROCx <- function(data_ROCx, metadata, plot_var, facet_var = NULL, plot_colors = NULL, xylim_one = FALSE){
  
  stopifnot(is.logical(xylim_one))
  stopifnot(length(data_ROCx) == nrow(metadata))
  stopifnot(plot_var %in% colnames(metadata))
  if(!is.null(facet_var))
    stopifnot(facet_var %in% colnames(metadata))
  
  metadata <- metadata[, c(plot_var, facet_var), drop = FALSE]
  
  if(!is.null(plot_colors)){
    stopifnot(nlevels(metadata[, plot_var]) == length(plot_colors))
  }else{
    plot_colors <- colorb(nlevels(metadata[, plot_var]))
  }

  
  ROClist <- lapply(1:length(data_ROCx), function(r){
    # r = 1
    
    ROC <- data_ROCx[[r]]$ROC
    
    ROC <- cbind(ROC, metadata[rep(r, nrow(ROC)), , drop = FALSE])
    
  })
  
  Xlist <- lapply(1:length(data_ROCx), function(r){
    # r = 1
    
    X <- data_ROCx[[r]]$X
    
    X <- cbind(X, metadata[rep(r, nrow(X)), , drop = FALSE])
    
  })
  
  X <- do.call(rbind, Xlist)
  ROC <- do.call(rbind, ROClist)

  ggp <- ggplot(data = ROC, aes_string(x = "FPR", y = "TPR", group = plot_var, colour = plot_var)) +
    theme_bw() +
    geom_line(size = 1.5, na.rm=TRUE) +
    theme(axis.text=element_text(size = 16), axis.title = element_text(size = 18, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 12), strip.text = element_text(size = 12)) +
    guides(colour = guide_legend(override.aes = list(size = 1.5, shape = NA), ncol = 3)) +
    geom_point(data = X, aes_string(x = "FPR", y = "TPR", group = plot_var, colour = plot_var), size = 10, shape = "X", na.rm=TRUE) 
  
  if(xylim_one)
    ggp <- ggp + coord_cartesian(xlim = c(-0.1, 1.1), ylim = c(-0.1, 1.1))
  
  ggp <- ggp + scale_color_manual(values = plot_colors)
  
  if(length(facet_var) == 1)
    ggp <- ggp + facet_wrap(reformulate(facet_var[1]))
  
  if(length(facet_var) == 2)
    ggp <- ggp + facet_grid(reformulate(facet_var[1], facet_var[2]))
  
  print(ggp)
  
}






