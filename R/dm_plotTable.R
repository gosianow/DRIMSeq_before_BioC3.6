

dm_plotTable <- function(pvalues, name = "pvalue", out_dir = NULL){
  
  df <- data.frame(pvalues = pvalues[!is.na(pvalues)])
  
  ggp <- ggplot(df, aes(x = pvalues)) +
    theme_bw() +
    xlab("p-values") +
    ylab("Frequency") +
    geom_histogram(binwidth = 0.01, fill = "deeppink4") +
    theme(axis.text = element_text(size=16), axis.title = element_text(size=18, face="bold"), plot.title = element_text(size=16, face="bold")) +
    coord_cartesian(xlim = c(-0.02, 1.02)) +
    geom_text(data = data.frame(x = Inf, y = Inf, label = paste0(nrow(df), " tests       ")), aes(x = x, y = y, label = label), hjust = 1, vjust = 3, size = 6)
  
  
  if(!is.null(out_dir))
  pdf(paste0(out_dir, "hist_", name, ".pdf"))
  
  print(ggp)
  
  
  # opar <- par() # make a copy of current settings
  # par(mar = c(5, 5, 4, 2) + 0.1, mgp = c(3, 1, 0)) # c(5, 4, 4, 2) + 0.1 # c(bottom, left, top, right)
  
  # hist(pvalues, breaks = 100, main = "", col = "hotpink", xlab = "p-values", cex.lab=1.5, cex.axis=1.5)
  # par(mar = opar$mar, mgp = opar$mgp) 
  
  if(!is.null(out_dir))
  dev.off()
  
  
}










