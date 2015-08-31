#' @import ggplot2

# counts = x@counts; out_dir = "~/"

dmDS_plotData <- function(counts, out_dir = NULL){

  tt <- width(counts)
  
  df <- data.frame(tt = tt)
  
  main <- paste0(length(tt), " genes \n", sum(tt) , " features ")
  
  ggp <- ggplot(df, aes(x = tt)) +
    theme_bw() +
    ggtitle(main) +
    xlab("Number of features per gene") +
    ylab("Frequency") +
    geom_histogram(binwidth = 1, fill = "forestgreen") +
    theme(axis.text = element_text(size=16), axis.title = element_text(size=18, face="bold"), plot.title = element_text(size=16, face="bold")) +
    coord_cartesian(xlim = c(0, max(tt) + 2)) 
  
  
  if(!is.null(out_dir))
  pdf(paste0(out_dir, "DM_hist_features.pdf"))
  
  print(ggp)
  
  # opar <- par()      # make a copy of current settings
  # par(mar = c(5, 5, 4, 2) + 0.1, mgp = c(3, 1, 0)) # c(5, 4, 4, 2) + 0.1 # c(bottom, left, top, right)
  
  # hist(tt, breaks = seq(0, max(tt), by = 1), col = "darkseagreen2", main = paste0(length(tt), " genes \n ", sum(tt) , " features "), xlab = "Number of features per gene", cex.lab=1.5, cex.axis=1.5)
  
  # par(mar = opar$mar, mgp = opar$mgp) 
  
  if(!is.null(out_dir))
  dev.off()
  
  
}










