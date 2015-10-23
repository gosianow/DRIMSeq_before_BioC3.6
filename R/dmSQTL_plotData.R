

dmSQTL_plotData <- function(counts, genotypes, blocks, out_dir = NULL){
  
  tt <- width(counts)
  
  df <- data.frame(tt = tt)
  binwidth <- ceiling(max(df$tt)/50)
  
  ggp1 <- ggplot(df, aes_string(x = "tt")) +
    theme_bw() +
    xlab("Number of features per gene") +
    ylab("Frequency") +
    geom_histogram(fill = "seagreen4", binwidth = binwidth) +
    theme(axis.text = element_text(size=16), axis.title = element_text(size=18, face="bold"), plot.title = element_text(size=18, face="bold")) +
    coord_cartesian(xlim = c(0, max(tt) + 2)) +
    geom_text(data = data.frame(x = Inf, y = Inf, label = paste0(length(tt), " genes   \n", sum(tt) , " features   ")), aes_string(x = "x", y = "y", label = "label"), hjust = 1, vjust = 2, size = 6)
  
  tt <- width(blocks)
  
  
  df <- data.frame(tt = tt)
  binwidth <- ceiling(max(df$tt)/50)
  
  ggp2 <- ggplot(df, aes_string(x = "tt")) +
    theme_bw() +
    xlab("Number of SNPs per gene") +
    ylab("Frequency") +
    geom_histogram(fill = "royalblue4", binwidth = binwidth) +
    theme(axis.text = element_text(size=16), axis.title = element_text(size=18, face="bold"), plot.title = element_text(size=18, face="bold")) +
    coord_cartesian(xlim = c(0, max(tt) + 2)) +
    geom_text(data = data.frame(x = Inf, y = Inf, label = paste0(length(tt), " genes   \n", sum(tt) , " SNPs   ")), aes_string(x = "x", y = "y", label = "label"), hjust = 1, vjust = 2, size = 6)
  
  
  tt <- width(genotypes)
  
  
  df <- data.frame(tt = tt)
  binwidth <- ceiling(max(df$tt)/50)
  
  ggp3 <- ggplot(df, aes_string(x = "tt")) +
    theme_bw() +
    xlab("Number of blocks per gene") +
    ylab("Frequency") +
    geom_histogram(fill = "mediumpurple4", binwidth = binwidth) +
    theme(axis.text = element_text(size=16), axis.title = element_text(size=18, face="bold"), plot.title = element_text(size=18, face="bold")) +
    coord_cartesian(xlim = c(0, max(tt) + 2)) +
    geom_text(data = data.frame(x = Inf, y = Inf, label = paste0(length(tt), " genes   \n", sum(tt) , " blocks   ")), aes_string(x = "x", y = "y", label = "label"), hjust = 1, vjust = 2, size = 6)
  
  
  if(!is.null(out_dir)){
    
    pdf(paste0(out_dir, "hist_features.pdf"))
    print(ggp1)
    dev.off()
    
    pdf(paste0(out_dir, "hist_snps.pdf"))
    print(ggp2)
    dev.off()
    
    pdf(paste0(out_dir, "hist_blocks.pdf"))
    print(ggp3)
    dev.off()
    
  }else{
    
    print(multiplot(ggp1, ggp2, ggp3, cols = 3))
    
  }
  
  
}










