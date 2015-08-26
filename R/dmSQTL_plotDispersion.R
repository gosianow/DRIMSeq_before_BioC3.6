#' @import ggplot2
dmSQTL_plotDispersion <- function(genewise_dispersion, mean_expression, nr_features, common_dispersion = numeric(), out_dir = NULL){
  
  w <- sapply(genewise_dispersion, length)
  
  mean_expression <- rep(mean_expression, w)
  nr_features <- rep(nr_features, w)
  
  df <- data.frame(mean_expression = log10(mean_expression + 1), dispersion = log10(unlist(genewise_dispersion)), nr_features = nr_features)
  
  df_quant <- quantile(na.omit(df$nr_features), probs = 0.95)
  
  
  ggp2 <- ggplot(df, aes(x = mean_expression, y = dispersion, colour = nr_features )) +
  theme_bw() +
  xlab("Log10(mean expression)") +
  ylab("Log10(dispersion)") +
  geom_point(size = 1.5, alpha = 0.5) +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=18, face="bold"), legend.title = element_text(size=16, face="bold"), legend.text = element_text(size = 14), legend.position = "top") +
  guides(colour = guide_colorbar(barwidth = 20, barheight = 0.5)) +
  scale_colour_gradient(limits = c(1, df_quant), breaks = seq(1, df_quant, 1), low = "dodgerblue2", high="firebrick2", name = "Number of features", na.value = "firebrick2")
  
  if(length(common_dispersion)){
    ggp2 <- ggp2 + geom_hline(yintercept = log10(common_dispersion), colour = "black", linetype = "dashed", size =  0.5)
  }

  if(!is.null(out_dir))
  pdf(paste0(out_dir, "dispersion_vs_mean.pdf"))
  
  print(ggp2)
  
  if(!is.null(out_dir))
  dev.off()
  
  
}










