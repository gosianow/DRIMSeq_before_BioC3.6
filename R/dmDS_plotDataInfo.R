# counts = x@counts; out_dir = "~/"; info = NULL

dmDS_plotDataInfo <- function(counts, out_dir = NULL, info = NULL){

  ### Plot the frequency of present features
  
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










