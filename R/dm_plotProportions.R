# counts must have rownames corresponding to features

dm_plotProportions <- function(counts, group, pi_full = NULL, pi_null = NULL, main = NULL, plot_type = "boxplot1", order = TRUE){
  
  labels <- labels_org <- factor(rownames(counts), levels = rownames(counts))
  group_counts <- table(group)
  
  prop_samp <- data.frame( feature_id = labels, prop.table(counts, 2), stringsAsFactors = FALSE) 
  
  if(!is.null(pi_full))
    prop_est_full <- data.frame(feature_id = labels, pi_full, stringsAsFactors = FALSE)
  if(!is.null(pi_null))
    prop_est_null <- data.frame(feature_id = labels, pi_null, stringsAsFactors = FALSE)
  
  #### order transcipts by decreasing proportion 
  if(order){
    labels <- labels[order(apply(aggregate(t(prop_samp[, -1]), by = list(group = group), median)[, -1], 2, max), decreasing = TRUE)]  
  }
  
  prop_samp <- melt(prop_samp, id.vars = "feature_id", variable.name = "sample_id", value.name = "proportion")  
  prop_samp$feature_id <- factor(prop_samp$feature_id, levels = labels)
  prop_samp$sample_id <- factor(prop_samp$sample_id)
  prop_samp$group <- rep(group, each = length(labels))
  
  if(!is.null(pi_full)){
    prop_est_full <- melt(prop_est_full, id.vars = "feature_id", variable.name = "group", value.name = "proportion")
    prop_est_full$feature_id <- factor(prop_est_full$feature_id, levels = labels)
    prop_est_full$group <- factor(rep(levels(group), each = length(labels)), levels = levels(group))
  }
  
  if(!is.null(pi_null)){
    colnames(prop_est_null) <- c("feature_id", "proportion")
    prop_est_null$feature_id <- factor(prop_est_null$feature_id, levels = labels)
    prop_est_null$group <- factor("null")
  }
  
  if(plot_type == "barplot"){
    
    values <- c(colorb(nlevels(group)), "grey65")
    names(values) <- c(levels(group), "null")
    order_prop_samp <- order(prop_samp$group, prop_samp$sample_id)
    prop_samp$sample_id <- factor(prop_samp$sample_id, levels = unique(prop_samp$sample_id[order_prop_samp]))
    width = 0.9
    
    if(!is.null(pi_full)){
      prop_est_full$xid <- as.numeric(prop_est_full$feature_id)
      prop_est_full$group_prop <- prop_est_full$group
      levels(prop_est_full$group_prop) <- group_counts/sum(group_counts)
      prop_est_full$group_prop <- as.numeric(as.character(prop_est_full$group_prop))
      prop_est_full$group_cumsum <- prop_est_full$group
      levels(prop_est_full$group_cumsum) <- cumsum(group_counts)/sum(group_counts)
      prop_est_full$group_cumsum <- as.numeric(as.character(prop_est_full$group_cumsum))
      prop_est_full$x <- prop_est_full$xid - width/2 + prop_est_full$group_cumsum * width -prop_est_full$group_prop/2
    }
    
    ggp <- ggplot() +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.text=element_text(size=14), axis.title = element_text(size=14, face="bold"), plot.title = element_text(size=14), legend.position = "right", legend.title = element_text(size = 14), legend.text = element_text(size = 14)) +
      ggtitle(main) +
      geom_bar(data = prop_samp, aes_string(x = "feature_id", y = "proportion", group = "sample_id", fill = "group"), stat = "identity", position = position_dodge(width = width)) +
      scale_fill_manual(name = "Groups", values = values, breaks = names(values)) +
      xlab("Features") +
      ylab("Proportions")
    
    if(!is.null(pi_null)){
      ggp <- ggp +
        geom_point(data = prop_est_null, aes_string(x = "feature_id", y = "proportion", fill = "group"), size = 3.5, shape = 22) 
    }
    
    if(!is.null(pi_full)){
      ggp <- ggp + 
        geom_point(data = prop_est_full, aes_string(x = "x", y = "proportion", fill = "group"), size = 3, shape = 23)
    }

    return(ggp)
    
  }
  
  if(plot_type == "boxplot1"){
    ### box plots with points
    
    values <- c(colorb(nlevels(group)), "grey65")
    names(values) <- c(levels(group), "null")
    
    ggp <- ggplot() +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.text=element_text(size=14), axis.title=element_text(size=14, face="bold"), plot.title = element_text(size=14), legend.position = "right", legend.title = element_text(size = 14), legend.text = element_text(size = 14)) +
      ggtitle(main) +     
      geom_jitter(data = prop_samp, aes_string(x = "feature_id", y = "proportion", fill = "group", colour = "group"), position = position_jitterdodge(dodge.width = 0.75), alpha = 0.5, size = 2, show_guide = FALSE, na.rm = TRUE) +
      geom_boxplot(data = prop_samp, aes_string(x = "feature_id", y = "proportion", colour = "group", fill = "group"), outlier.size = 0, alpha = 0.2, lwd = 0.5) +
      scale_fill_manual(name = "Groups", values = values, breaks = names(values)) +
      scale_colour_manual(name = "Groups", values = values, breaks = names(values)) +
      xlab("Features") +
      ylab("Proportions")
    
    if(!is.null(pi_null)){
      ggp <- ggp +
        geom_point(data = prop_est_null, aes_string(x = "feature_id", y = "proportion", fill = "group"), size = 3.5, shape = 22) +
        guides(colour=FALSE)
    }
    
    if(!is.null(pi_full)){
      ggp <- ggp + 
        geom_point(data = prop_est_full, aes_string(x = "feature_id", y = "proportion", fill = "group"), position = position_jitterdodge(jitter.width = 0, jitter.height = 0), size = 3, shape = 23) +
        guides(colour=FALSE)
    }

    return(ggp)
  
  }
  
  if(plot_type == "boxplot2"){
    ### box plots per group
    values <- colorb(length(labels_org))
    names(values) <- labels_org
    
    ggp <- ggplot() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 0, vjust = 0.5), axis.text=element_text(size=14), axis.title=element_text(size=14, face="bold"), plot.title = element_text(size=14), panel.grid.major = element_blank(), legend.title = element_text(size = 14), legend.text = element_text(size = 14)) +
      geom_vline(xintercept = seq(1, nlevels(group) - 1, 1) + 0.5, color = "gray90") +
      ggtitle(main) +     
      geom_boxplot(data = prop_samp, aes_string(x = "group", y = "proportion", fill = "feature_id"), width = 1) + 
      scale_fill_manual(name = "Features", values = values) +
      scale_x_discrete(labels = paste0(names(group_counts), " (", group_counts, ")" ), name="") +
      guides(fill = guide_legend(nrow = 20)) +
      xlab("Groups") +
      ylab("Proportions")
    
    if(!is.null(pi_full)){
      ggp <- ggp + 
        geom_point(data = prop_est_full, aes_string(x = "group", y = "proportion", fill = "feature_id"), position = position_jitterdodge(jitter.width = 0, jitter.height = 0), size = 3, shape = 23, colour = "black")
    }
    
    return(ggp)
    
  }
  
  if(plot_type == "lineplot"){
    ### line plots
    values <- c(colorb(nlevels(group)), "grey65")
    names(values) <- c(levels(group), "null")
    
    ggp <- ggplot() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.text=element_text(size=14), axis.title =element_text(size=14, face="bold"), plot.title = element_text(size=14), legend.position = "right", legend.title = element_text(size = 14), legend.text = element_text(size = 14)) +
      ggtitle(main) +
      geom_line(data = prop_samp, aes_string(x = "feature_id", y = "proportion", group = "sample_id", colour = "group"), size = 1.1) +
      scale_fill_manual(name = "Groups", values = values, breaks = names(values)) +
      scale_colour_manual(name = "Groups", values = values, breaks = names(values)) +
      xlab("Features") +
      ylab("Proportions")
    
    if(!is.null(pi_null)){
      ggp <- ggp +
        geom_point(data = prop_est_null, aes_string(x = "feature_id", y = "proportion", fill = "group"), size = 3.5, shape = 22, show_guide = FALSE) +
        guides(colour=FALSE)
    }
    
    if(!is.null(pi_full)){
      ggp <- ggp + 
        geom_point(data = prop_est_full, aes_string(x = "feature_id", y = "proportion", group = "group", fill = "group"), size = 3, shape = 23) +
        guides(colour=FALSE)
    }
  
    return(ggp)

  }
  
  if(plot_type == "ribbonplot" & !is.null(pi_full)){
    
    values <- colorb(length(labels_org))
    names(values) <- labels_org
    breaks <- labels
    width  <- 0.5
    prop_est_full_order <- prop_est_full[order(prop_est_full$group), ]
    
    if(order == TRUE){
      prop_est_full_order <- prop_est_full[order(prop_est_full$group, prop_est_full$proportion), ]
      breaks = rev(prop_est_full_order[prop_est_full_order$group == levels(prop_est_full_order$group)[1], "feature_id"])
    }

    ### get ribbons!!!
    gr <- list()
    
    for (i in 1:(nlevels(group) - 1)){
      # i = 2
      prop_est_full_ribbon <- prop_est_full_order[prop_est_full_order$group %in% levels(prop_est_full_order$group )[c(i, i+1)], ]
      prop_est_full_ribbon$group <- factor(prop_est_full_ribbon$group)
      prop_est_full_ribbon$cumsum <- matrix(t(aggregate(prop_est_full_ribbon[,"proportion"], by = list(group = prop_est_full_ribbon$group), cumsum)[, -1]), ncol = 1)
      prop_est_full_ribbon$offset <- c(width/2, -width/2)[as.numeric(prop_est_full_ribbon$group)]
      prop_est_full_ribbon$xid <- i - 1
      prop_est_full_ribbon$x <- as.numeric(prop_est_full_ribbon$group) + prop_est_full_ribbon$offset + prop_est_full_ribbon$xid
      prop_est_full_ribbon$ymin <- prop_est_full_ribbon$cumsum - prop_est_full_ribbon$proportion
      prop_est_full_ribbon$ymax <- prop_est_full_ribbon$cumsum
      
      gr[[i]] <- geom_ribbon(data = prop_est_full_ribbon, aes_string(x = "x", ymin = "ymin", ymax = "ymax", group = "feature_id", fill = "feature_id"), alpha = 0.3) 
      
    }

    ggp <- ggplot() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 0, vjust = 0.5), axis.text=element_text(size=14), axis.title=element_text(size=14, face="bold"), plot.title = element_text(size=14), legend.title = element_text(size = 14), legend.text = element_text(size = 14)) +
      ggtitle(main) +    
      coord_cartesian(ylim = c(-0.1, 1.1)) + 
      coord_cartesian(ylim = c(-0.1, 1.1)) +
      scale_fill_manual(name = "Features", values = values, breaks = breaks) +
      scale_x_discrete(labels = paste0(names(group_counts), " (", group_counts, ")" ), name="") +
      guides(fill = guide_legend(nrow = 25)) +
      xlab("Groups") +
      ylab("Estimated proportions") +
      geom_bar(data = prop_est_full_order, aes_string(x = "group", y = "proportion", fill = "feature_id"), stat = "identity", width = width, position="stack") + gr
    
    return(ggp)
    
  }
  
}












