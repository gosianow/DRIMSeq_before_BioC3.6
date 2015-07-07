colorb <- function(n) {
  colorRampPalette(c("#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))(n)
}

colorl <- function(n){
  colorRampPalette(c("#CDC0B0", "#66CDAA", "#838B8B", "#0000EE", "#CD3333", "#7AC5CD", "#66CD00", "#CD661D", "#6495ED", "#FFB90F", "#9A32CD", "#EE6A50", "#CAFF70", "#E9967A", "#483D8B", "#CD1076", "#00CED1", "#EE2C2C", "#1C86EE", "#EEC900", "#00CD00", "#EE6AA7", "#ADFF2F", "#0000CD", "#8470FF", "#F08080", "#3CB371", "#FF34B3", "#9ACD32", "#FFA07A", "#CD3700", "#7EC0EE", "#8B668B", "#43CD80", "#7EC0EE", "#7D26CD", "#EEEE00", "#EE5C42", "#00CD66", "#6C7B8B", "#EE82EE", "#EED8AE", "#9ACD32", "#00C5CD", "#CD3278"))(n)
}

dm_plotProportion <- function(counts, group, sample_id, fit_full = NULL, fit_null = NULL, main = NULL, plot_type = "boxplot1", order = TRUE){

  labels <- rownames(counts)
  group_counts <- table(group)
  
  prop_samp <- data.frame( feature_id = labels, prop.table(counts, 2), stringsAsFactors = FALSE) 
  
  if(!is.null(fit_full))
  prop_est_full <- data.frame(feature_id = labels, fit_full$pi, stringsAsFactors = FALSE)
  if(!is.null(fit_null))
  prop_est_null <- data.frame(feature_id = labels, fit_null$pi, stringsAsFactors = FALSE)
  
  #### order transcipts by decreasing proportion 
  if(order){
    labels <- labels[order(apply(aggregate(t(prop_samp[, -1]), by = list(group = group), median)[, -1], 2, max), decreasing = TRUE)]  
  }
  
  
  prop_samp <- melt(prop_samp, id.vars = "feature_id", variable.name = "sample_id", value.name = "proportion")  
  prop_samp$feature_id <- factor(prop_samp$feature_id, levels = labels)
  prop_samp$sample_id <- factor(prop_samp$sample_id)
  prop_samp$group <- rep(group, each = length(labels))
  
  if(!is.null(fit_full)){
    prop_est_full <- melt(prop_est_full, id.vars = "feature_id", variable.name = "group", value.name = "proportion")
    prop_est_full$feature_id <- factor(prop_est_full$feature_id, levels = labels)
    prop_est_full$group <- prop_est_full$group
  }

  if(!is.null(fit_null)){
    colnames(prop_est_null) <- c("feature_id", "proportion")
    prop_est_null$feature_id <- factor(prop_est_null$feature_id, levels = labels)
    prop_est_null$group <- factor("null")
  }

  
  
  if(plot_type == "barplot"){

    values <- c(colorb(nlevels(group)), "#E69F00")
    names(values) <- c(levels(group), "null")
    
    ggp <- ggplot() +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.text=element_text(size=12), axis.title = element_text(size=12, face="bold"), plot.title = element_text(size=10)) +
    ggtitle(main) +
    geom_bar(data = prop_samp, aes(x = feature_id, y = proportion, group = sample_id, fill = group), stat = "identity", position = position_dodge()) +
    scale_fill_manual(name = "Groups", values = values, breaks = names(values)) +
    xlab("Features") +
    ylab("Proportions") +
    coord_cartesian(ylim = c(-0.1, 1.1))
    
    if(!is.null(fit_full)){
      ggp <- ggp + 
      geom_point(data = prop_est_full, aes(x = feature_id, y = proportion, fill = group), position = position_jitterdodge(jitter.width = 0, jitter.height = 0 ), size = 3, shape = 23)
      
    }
    
    if(!is.null(fit_null)){
      ggp <- ggp +
      geom_point(data = prop_est_null, aes(x = feature_id, y = proportion, fill = group), size = 3, shape = 23) 
    }
    
    
    return(ggp)
    
  }
  
  if(plot_type == "boxplot1"){
    ### box plots with points
    
    values <- c(colorb(nlevels(group)), "#E69F00")
    names(values) <- c(levels(group), "null")
    
    ggp <- ggplot() +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.text=element_text(size=12), axis.title=element_text(size=12, face="bold"), plot.title = element_text(size=10)) +
    ggtitle(main) +     
    geom_jitter(data = prop_samp, aes(x = feature_id, y = proportion, fill = group, colour = group), position = position_jitterdodge(dodge.width = 0.75), alpha = 0.5, size = 2, show_guide = FALSE) +
    geom_boxplot(data = prop_samp, aes(x = feature_id, y = proportion, colour = group), fill = "white", outlier.size = NA, alpha = 0, lwd = 0.5) +
    coord_cartesian(ylim = c(-0.1, 1.1))  +
    scale_fill_manual(name = "Groups", values = values, breaks = names(values)) +
    scale_colour_manual(name = "Groups", values = values, breaks = names(values)) +
    xlab("Features") +
    ylab("Proportions")
    
    if(!is.null(fit_full)){
      ggp <- ggp + 
      geom_point(data = prop_est_full, aes(x = feature_id, y = proportion, fill = group), position = position_jitterdodge(jitter.width = 0, jitter.height = 0), size = 3, shape = 23) +
      guides(colour=FALSE)
      
    }
    
    if(!is.null(fit_null)){
      ggp <- ggp +
      geom_point(data = prop_est_null, aes(x = feature_id, y = proportion, fill = group), size = 3, shape = 23) +
      guides(colour=FALSE)
    }

    
    return(ggp)
    
    
  }
  
  
  if(plot_type == "boxplot2"){

    ### box plots per group
    ggp <- ggplot() +
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5), axis.text=element_text(size=12), axis.title=element_text(size=12, face="bold"), plot.title = element_text(size=10), panel.grid.major = element_blank()) +
    geom_vline(xintercept = seq(1, nlevels(group) - 1, 1) + 0.5, color="white") +
    ggtitle(main) +     
    geom_boxplot(data = prop_samp, aes(x = group, y = proportion, fill = feature_id), width = 1) + 
    coord_cartesian(ylim = c(-0.1, 1.1)) +
    scale_fill_manual(name = "Features", values = colorl(length(labels))) +
    scale_x_discrete(labels = paste0(names(group_counts), " (", group_counts, ")" ), name="") +
    guides(fill = guide_legend(nrow = 25)) +
    xlab("Groups") +
    ylab("Proportions")
    
    
    if(!is.null(fit_full)){
      ggp <- ggp + 
      geom_point(data = prop_est_full, aes(x = group, y = proportion, fill = feature_id), position = position_jitterdodge(jitter.width = 0, jitter.height = 0), size = 3, shape = 23, colour = "black")
      
    }
    

    return(ggp)
    
  }
  
  if(plot_type == "lineplot"){
    ### line plots
    
    values <- c(colorb(nlevels(group)), "#E69F00")
    names(values) <- c(levels(group), "null")
    
    ggp <- ggplot() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.text=element_text(size=12), axis.title =element_text(size=12, face="bold"), plot.title = element_text(size=10)) +
    ggtitle(main) +
    geom_line(data = prop_samp, aes(x = feature_id, y = proportion, group = sample_id, colour = group), size = 1.1) +
    coord_cartesian(ylim = c(-0.1, 1.1)) +
    scale_fill_manual(name = "Groups", values = values, breaks = names(values)) +
    scale_colour_manual(name = "Groups", values = values, breaks = names(values)) +
    xlab("Features") +
    ylab("Proportions")
    
    
    if(!is.null(fit_full)){
      ggp <- ggp + 
      geom_point(data = prop_est_full, aes(x = feature_id, y = proportion, group = group, fill = group ), size = 3, shape = 23) +
      guides(colour=FALSE)
      
    }
    
    if(!is.null(fit_null)){
      ggp <- ggp +
      geom_point(data = prop_est_null, aes(x = feature_id, y = proportion, fill = group), size = 3, shape = 23, show_guide = FALSE) +
      guides(colour=FALSE)
    }
    
    
    return(ggp)
    
    
  }
  
  if(plot_type == "ribbonplot" & !is.null(fit_full)){

### Work only for two groups!!!

width  <- 0.5
prop_est_full_order <- prop_est_full[order(prop_est_full$group, prop_est_full$proportion), ]

prop_est_full_ribbon <- prop_est_full_order
prop_est_full_ribbon$cumsum <- matrix(t(aggregate(prop_est_full_order[,"proportion"], by = list(group = prop_est_full_order$group), cumsum)[, -1]), ncol = 1)
prop_est_full_ribbon$offset <- c(width/2, -width/2)[as.numeric(prop_est_full_ribbon$group)]


ggp <- ggplot() +
theme_bw() +
theme(axis.text.x = element_text(angle = 0, vjust = 0.5), axis.text=element_text(size=12), axis.title=element_text(size=12, face="bold"), plot.title = element_text(size=10)) +
ggtitle(main) +    
coord_cartesian(ylim = c(-0.1, 1.1)) + 
coord_cartesian(ylim = c(-0.1, 1.1)) +
scale_fill_manual(name = "Features", values = colorl(length(labels))) +
scale_x_discrete(labels = paste0(names(group_counts), " (", group_counts, ")" ), name="") +
guides(fill = guide_legend(nrow = 25)) +
xlab("Groups") +
ylab("Proportions") +
geom_bar(data = prop_est_full_order, aes(x = group, y = proportion, fill = feature_id), stat = "identity", width = width, position="stack") +
geom_ribbon(data = prop_est_full_ribbon, aes(x = as.numeric(group) + offset, ymin = cumsum - proportion, ymax = cumsum, group = feature_id, fill = feature_id), alpha = 0.7) 


    # pdf(paste0(out_dir, "proportions_", gsub(pattern = "\\.", replacement = "_" , gene), ".pdf"))
    # print(ggp)
    # dev.off()

    return(ggp)

  }






}












