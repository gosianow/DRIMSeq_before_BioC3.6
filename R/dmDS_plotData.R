

dmDS_plotData <- function(counts, out_dir = NULL, info = NULL){

  tt <- width(counts@partitioning)
  
  if(!is.null(out_dir))
  pdf(paste0(out_dir, "hist_feature_number.pdf"))
  
  opar <- par()      # make a copy of current settings
  par(mar = c(5, 5, 4, 2) + 0.1, mgp = c(3, 1, 0)) # c(5, 4, 4, 2) + 0.1 # c(bottom, left, top, right)
  
  hist(tt, breaks = seq(0, max(tt), by = 1), col = "darkseagreen2", main = paste0(length(tt), " genes \n ", sum(tt) , " features "), xlab = "Number of features per gene", cex.lab=1.5, cex.axis=1.5)
  
  par(mar = opar$mar, mgp = opar$mgp) 
  if(!is.null(out_dir))
  dev.off()
  
  if(!is.null(info)){

    info_spl <- split(info, info$gene_id)
    
    genes <- names(counts)
    
    tas <- table(unlist(lapply(names(info_spl), function(g){
      # g = "FBgn0004636"
      feature <- info_spl[[g]][, "feature_id"]
      
      if(! g %in% genes)
        return(NA)
      
      feature_in <- sum(rownames(counts[[g]]) %in% feature)
      
      return(feature_in)
      
    })), useNA = "always")
    
    
    tas <- tas[c(length(tas), 1:(length(tas) -1))]
    names(tas)[1] <- "NG"
    
    colors <- c("darkred", "darkred", rep("grey", 200))
    names(colors) <- c("NG", "0", 1:200)
    colors <- colors[names(tas)]
    
    if(!is.null(out_dir))
    pdf(paste0(out_dir, "hist_info_filtering.pdf"), width = 7, height = 7)
    opar <- par() 
    par(mar = c(5, 5, 4, 2) + 0.1, mgp = c(3, 1, 0)) # c(5, 4, 4, 2) + 0.1 # c(bottom, left, top, right)
    xx <- barplot(tas, xlab = "Number of DS features left within DS gene", ylab = "Number of DS genes", col = colors, cex.lab=1.5)
    text(x = xx, y = as.numeric(tas), label = as.numeric(tas), pos = 3)
    par(mar = opar$mar, mgp = opar$mgp) 
    if(!is.null(out_dir))
    dev.off()
    
  }
  
  
}










