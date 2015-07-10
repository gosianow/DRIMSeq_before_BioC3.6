

dmSQTL_plotData <- function(data, out_dir = "./"){
  

  tt <- sapply(data@counts, nrow)
  
  pdf(paste0(out_dir, "hist_feature_number.pdf"))
  
  #opar <- par()      # make a copy of current settings
  par(mar = c(5, 5, 4, 2) + 0.1, mgp = c(3, 1, 0)) # c(5, 4, 4, 2) + 0.1 # c(bottom, left, top, right)
  
  hist(tt, breaks = seq(0, max(tt), by = 1), col = "darkseagreen2", main = paste0(length(tt), " genes \n ", sum(tt) , " features "), xlab = "Number of features per gene", cex.lab=1.5, cex.axis=1.5, cex.main = 1.5)
  
  # par(mar = opar$mar, mgp = opar$mgp) 
  
  dev.off()
  


  tt <- sapply(data@genotypes, nrow)
  
  pdf(paste0(out_dir, "hist_snp_number.pdf"))
  
  #opar <- par()      # make a copy of current settings
  par(mar = c(5, 5, 4, 2) + 0.1, mgp = c(3, 1, 0)) # c(5, 4, 4, 2) + 0.1 # c(bottom, left, top, right)
  
  hist(tt, breaks = 100, col = "darkturquoise", main = paste0(length(tt), " genes \n ", sum(tt) , " SNPs "), xlab = "Number of SNPs per gene", cex.lab=1.5, cex.axis=1.5, cex.main = 1.5)
  
  # par(mar = opar$mar, mgp = opar$mgp) 
  
  dev.off()
  
  
  
}










