#' TPR versus achieved FDR plot
#' 
#' TPR versus achieved FDR plot for a list of result data frames.
#' 
#' @inheritParams plotROCx
#' @param thresholds Thresholds for the FDR.
#' @return 
#' \code{calculateTPRFDR} returns a list of data frames with FDR thresholds and the corresponding TPRs and achieved FDRs.
#' 
#' \code{plotTPRFDR} returns a plot with TPR versus achieved FDR curves. TPR and achieved FDR are marked for each threshold with a circle. If the achieved FDR is lower than the corresponding threshold, the circle is filled with colour, otherwise, it is white inside.
#' 
#'  @author Malgorzata Nowicka
#'  @name plotTPRFDR
#' @export
calculateTPRFDR <- function(results, status, thresholds = c(0.01, 0.05, 0.1)){
  
  stopifnot(all(c("gene_id", "status") %in% colnames(status)))
  
  status <- status[complete.cases(status[, c("gene_id", "status")]), , drop = FALSE]
  
  TPRFDRlist <- list()
  
  for(m in 1:length(results)){
    # m = 1
    # print(m)
    
    if(! "adj_pvalue" %in% colnames(results[[m]])){
      
      message("No 'adj_pvalue' column in results ", m)
      TPRFDRlist[[m]] <- data.frame(threshold = NA, FDR = NA, TPR = NA)
      
      
    }else{
      
      sc <- merge(status, results[[m]], by = "gene_id", all.x = TRUE)
      
      # mm <- match(status$gene_id, results[[m]]$gene_id)
      
      apvs <- sc$adj_pvalue
      apvs[is.na(apvs)] <- 1
      sts <- sc$status
      
      q <- length(thresholds)
      TPR <- rep(0, q)
      FDR <- rep(0, q)
      
      for(i in 1:q){
        # i=1
        sts_est <- as.numeric(apvs < thresholds[i])
        
        TP <- sum(sts==1 & sts_est==1)
        FP <- sum(sts==0 & sts_est==1)
        FN <- sum(sts==1 & sts_est==0)
        
        TPR[i] <- TP/(TP+FN)
        FDR[i] <- FP/(FP+TP)
        
      }  
      
      TPRFDRlist[[m]] <- data.frame(threshold = thresholds, FDR = FDR, TPR = TPR)
      
    }
    
    
  }
  
  return(TPRFDRlist)
  
  
}



################################################################################

#' @param data_TPRFDR List structured like the output of \code{calculateTPRFDR}.
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
#' data_TPRFDR <- calculateTPRFDR(results, status, 
#'    thresholds = c(0.01, 0.05, 0.1))
#' 
#' plotTPRFDR(data_TPRFDR, metadata, plot_var = "method", facet_var = NULL, 
#'    plot_colors = NULL, xylim_one = TRUE)
#' 
#' @seealso \code{\link{dataDS_status}}, \code{\link{dataDS_dmDStest}}, \code{\link{plotROCx}}, \code{\link{plotVenn}}
#' @rdname plotTPRFDR
#' @export
plotTPRFDR <- function(data_TPRFDR, metadata, plot_var, facet_var = NULL, plot_colors = NULL, xylim_one = FALSE){
  
  stopifnot(is.logical(xylim_one))
  stopifnot(length(data_TPRFDR) == nrow(metadata))
  stopifnot(plot_var %in% colnames(metadata))
  if(!is.null(facet_var))
    stopifnot(facet_var %in% colnames(metadata))
  
  metadata <- metadata[, c(plot_var, facet_var), drop = FALSE]
  
  if(!is.null(plot_colors)){
    stopifnot(nlevels(metadata[, plot_var]) == length(plot_colors))
  }else{
    plot_colors <- colorb(nlevels(metadata[, plot_var]))
  }
  
  TPRFDRlist <- lapply(1:nrow(metadata), function(r){
    # r = 1
    
    TPRFDR <- data_TPRFDR[[r]]
    
    TPRFDR <- cbind(TPRFDR, metadata[rep(r, nrow(TPRFDR)), , drop = FALSE])
    
  })

  TPRFDR <- do.call(rbind, TPRFDRlist)
  
  TPRFDR$white <- ifelse(TPRFDR$FDR <= TPRFDR$threshold, NA, TPRFDR$TPR)
  
  pointsize <- 2.5
  
  ggp <- ggplot(data = TPRFDR, aes_string(x = "FDR", y = "TPR", group = plot_var, colour = plot_var)) +
    theme_bw() +
    xlab("FDR") +
    geom_line(size = 1.5, na.rm=TRUE) +
    geom_vline(aes_string(xintercept = "threshold"), linetype = "dashed") + 
    geom_point(size = pointsize + 1, shape = 19, na.rm=TRUE) + 
    geom_point(aes_string(y = "white"), size = pointsize, shape = 21, fill = "white", na.rm=TRUE) + 
    theme(axis.text=element_text(size = 16), axis.title = element_text(size = 18, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 12), strip.text = element_text(size = 12)) +
    guides(colour = guide_legend(override.aes = list(size = 1.5, shape = NA), ncol = 3)) 
  
  if(xylim_one)
    ggp <- ggp + coord_cartesian(xlim = c(-0.1, 1.1), ylim = c(-0.1, 1.1))
  
    ggp <- ggp + scale_color_manual(values = plot_colors)
  
  if(length(facet_var) == 1)
    ggp <- ggp + facet_wrap(reformulate(facet_var[1]))
  
  if(length(facet_var) == 2)
    ggp <- ggp + facet_grid(reformulate(facet_var[1], facet_var[2]))
  
  print(ggp)
  
}



