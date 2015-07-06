

dmDS_plotTable <- function(table, out_dir = "./"){
  

  tt <- sapply(data@counts, nrow)
  
  pdf(paste0(out_dir, "hist_pvalue.pdf"))
  
  #opar <- par()      # make a copy of current settings
  par(mar = c(5, 5, 4, 2) + 0.1, mgp = c(3, 1, 0)) # c(5, 4, 4, 2) + 0.1 # c(bottom, left, top, right)
  
  hist(table$pvalue, breaks = 100, main = "", col = "hotpink", xlab = "p-values", cex.lab=1.5, cex.axis=1.5)

  # par(mar = opar$mar, mgp = opar$mgp) 
  
  dev.off()
  
}










