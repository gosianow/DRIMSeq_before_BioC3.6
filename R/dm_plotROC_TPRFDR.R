

#' Calculate data for ROC plot
#' 
#' Calculate TPR and FPR for ROC plot
#' 
#' @param results list of data.frames with "gene_id", "pvalue" and "adj_pvalue" columns
#' @param status data.frame with "gene_id" and "status" columns. "status" must consists of 0 and 1 (or FALSE, TRUE)
#' @export
calculate_ROCx <- function(results, status){

  status <- status[complete.cases(status[, c("gene_id", "status")]), , drop = FALSE]
  
  P <- sum(status$status == 1)
  N <- sum(status$status == 0)
  
  Xlist <- list()
  ROClist <- list()
  
  for(m in 1:length(results)){
    # m = 2
    print(m)
    
    if(! "pvalue" %in% colnames(results[[m]])){
      
      message("No 'pvalue' column in results ", m)
      ROClist[[m]] <- data.frame(FPR = NA, TPR = NA)
      Xlist[[m]] <- data.frame(FPR = NA, TPR = NA) 
      
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
      ROClist[[m]] <- data.frame(FPR = FPRv, TPR = TPRv)
      
      ### what is the TPR if the error was controlled 
      TPR <- sum(apvs[sts == 1] < 0.05) / P    
      Xlist[[m]] <- data.frame(FPR = approx(TPRv, FPRv, xout=TPR)$y, TPR = TPR)    
      
      
    }
    

  }
  
  return(list(ROClist = ROClist, Xlist = Xlist))
   
}


################################################################################

#' @export
plot_ROCx <- function(data_ROCx, split_levels, plot_levels, facet_levels, plot_colors = NULL, xylim_one = FALSE){
  
  split_levels <- split_levels[, c(plot_levels, facet_levels), drop = FALSE]
  
  ROClist <- lapply(1:nrow(split_levels), function(r){
    # r = 1
    
    ROC <- data_ROCx$ROClist[[r]]
    
    ROC <- cbind(ROC, split_levels[rep(r, nrow(ROC)), , drop = FALSE])
    
  })
  
  Xlist <- lapply(1:nrow(split_levels), function(r){
    # r = 1
    
    X <- data_ROCx$Xlist[[r]]
    
    X <- cbind(X, split_levels[rep(r, nrow(X)), , drop = FALSE])
    
  })
  
  X <- do.call(rbind, Xlist)
  ROC <- do.call(rbind, ROClist)
    
    
  ggp <- ggplot(data = ROC, aes_string(x = "FPR", y = "TPR", group = plot_levels, colour = plot_levels)) +
  theme_bw() +
  geom_line(size = 1.5, na.rm=TRUE) +
  theme(axis.text=element_text(size = 16), axis.title = element_text(size = 18, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 12), strip.text = element_text(size = 12)) +
  guides(colour = guide_legend(override.aes = list(size = 1.5, shape = NA), ncol = 3)) +
  geom_point(data = X, aes_string(x = "FPR", y = "TPR", group = plot_levels, colour = plot_levels), size = 8, shape = "X", na.rm=TRUE) 
  
  if(xylim_one)
  ggp <- ggp + coord_cartesian(xlim = c(-0.1, 1.1), ylim = c(-0.1, 1.1))
      
  if(!is.null(plot_colors) && nlevels(ROC[, plot_levels]) == length(plot_colors))     
    ggp <- ggp + scale_color_manual(values = plot_colors)
  else
    ggp <- ggp + scale_color_manual(values = colorb(nlevels(TPRFDR[, plot_levels])))
      
  if(length(facet_levels) == 1)
    ggp <- ggp + facet_wrap(reformulate(facet_levels[1]))
      
  if(length(facet_levels) == 2)
    ggp <- ggp + facet_grid(reformulate(facet_levels[1], facet_levels[2]))

    
  # pdf("./ROC.pdf")
   
  print(ggp)

  # dev.off()
  
  
}



################################################################################

#' Calculate data for TPR versus achieved FDR plot
#' 
#' Calculate TPR and FDR.
#' 
#' @param results list of data.frames with "gene_id" and "adj_pvalue" columns
#' @param status data.frame with "gene_id" and "status" columns. "status" must consists of 0 and 1 (or FALSE, TRUE)
#' @export
calculate_TPRFDR <- function(results, status, thresholds = c(0.01, 0.05, 0.1)){

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

#' @export
plot_TPRFDR <- function(data_TPRFDR, split_levels, plot_levels, facet_levels, plot_colors = NULL, xylim_one = FALSE){
  
  split_levels <- split_levels[, c(plot_levels, facet_levels), drop = FALSE]
  
  
  TPRFDRlist <- lapply(1:nrow(split_levels), function(r){
    # r = 1
    
    TPRFDR <- data_TPRFDR[[r]]
    
    TPRFDR <- cbind(TPRFDR, split_levels[rep(r, nrow(TPRFDR)), , drop = FALSE])
    
  })
  

  TPRFDR <- do.call(rbind, TPRFDRlist)
  
  TPRFDR$white <- ifelse(TPRFDR$FDR <= TPRFDR$threshold, NA, TPRFDR$TPR)
  
  pointsize <- 2.5

  ggp <- ggplot(data = TPRFDR, aes_string(x = "FDR", y = "TPR", group = plot_levels, colour = plot_levels)) +
  theme_bw() +
  xlab("FDR") +
  geom_line(size = 1.5, na.rm=TRUE) +
  geom_vline(aes(xintercept = threshold), linetype = "dashed") + 
  geom_point(size = pointsize + 1, shape = 19, na.rm=TRUE) + 
  geom_point(aes_string(y = "white"), size = pointsize, shape = 21, fill = "white", na.rm=TRUE) + 
  theme(axis.text=element_text(size = 16), axis.title = element_text(size = 18, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 12), strip.text = element_text(size = 12)) +
  guides(colour = guide_legend(override.aes = list(size = 1.5, shape = NA), ncol = 3)) 
  
  if(xylim_one)
  ggp <- ggp + coord_cartesian(xlim = c(-0.1, 1.1), ylim = c(-0.1, 1.1))
      
  if(!is.null(plot_colors) && nlevels(TPRFDR[, plot_levels]) == length(plot_colors))   
    ggp <- ggp + scale_color_manual(values = plot_colors)
  else
    ggp <- ggp + scale_color_manual(values = colorb(nlevels(TPRFDR[, plot_levels])))
      
  if(length(facet_levels) == 1)
    ggp <- ggp + facet_wrap(reformulate(facet_levels[1]))
      
  if(length(facet_levels) == 2)
    ggp <- ggp + facet_grid(reformulate(facet_levels[1], facet_levels[2]))

    
  # pdf("./TPRFDR.pdf")
   
  print(ggp)

  # dev.off()
  

}



