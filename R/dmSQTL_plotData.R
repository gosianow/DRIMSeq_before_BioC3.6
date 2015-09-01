

dmSQTL_plotData <- function(counts, genotypes, out_dir = NULL){
  
  tt <- width(counts)
  
  df <- data.frame(tt = tt)
  
  ggp <- ggplot(df, aes(x = tt)) +
    theme_bw() +
    ggtitle(paste0(length(tt), " genes \n ", sum(tt) , " features ")) +
    xlab("Number of features per gene") +
    ylab("Frequency") +
    geom_histogram(binwidth = 1, fill = "seagreen4") +
    theme(axis.text = element_text(size=16), axis.title = element_text(size=18, face="bold"), plot.title = element_text(size=16, face="bold")) +
    coord_cartesian(xlim = c(0, max(tt) + 2))
  
  
  if(!is.null(out_dir))
  pdf(paste0(out_dir, "DM_hist_features.pdf"))
  
  print(ggp)
  
  if(!is.null(out_dir))
  dev.off()




  tt <- width(genotypes)
  
  
  df <- data.frame(tt = tt)
  
  ggp <- ggplot(df, aes(x = tt)) +
    theme_bw() +
    ggtitle(paste0(length(tt), " genes \n ", sum(tt) , " SNPs ")) +
    xlab("Number of SNPs per gene") +
    ylab("Frequency") +
    geom_histogram(binwidth = 10, fill = "royalblue4") +
    theme(axis.text = element_text(size=16), axis.title = element_text(size=18, face="bold"), plot.title = element_text(size=16, face="bold")) +
    coord_cartesian(xlim = c(0, max(tt) + 2))
  
  
  if(!is.null(out_dir))
  pdf(paste0(out_dir, "DM_hist_snps.pdf"))
  
  print(ggp)
  
  if(!is.null(out_dir))
  dev.off()
  
  
  
  
}










