##############################################################################
# multiple group ffting 
##############################################################################
# model = c("full", "null")[1]; prop_mode=c("constrOptim2", "constrOptim2G", "FisherScoring")[2]; prop_tol = 1e-12; verbose = FALSE; BPPARAM = MulticoreParam(workers=10)

dmDS_fitOneModel <- function(counts, samples, dispersion, model = c("full", "null")[1], prop_mode=c("constrOptim2", "constrOptim2G", "FisherScoring")[2], prop_tol = 1e-12, verbose = FALSE, BPPARAM = MulticoreParam(workers=1)){
  
  gene_list <-  names(counts)
	
  if(length(dispersion) == 1){
    gamma0 <- rep(dispersion, length(gene_list))
		names(gamma0) <- gene_list
  } else {
    gamma0 <- dispersion
  }
  
  
  switch(model, 
         
         full = {
           cat("Fitting full model.. \n")
           group <- samples$group
           ngroups <- nlevels(group)
           lgroups <- levels(group)
           
           igroups <- lapply(lgroups, function(gr){which(group == gr)})
           names(igroups) <- lgroups
           
           
           time <- system.time(ff <- BiocParallel::bplapply(gene_list, function(g){  
             # g = "FBgn0000008"
             # cat("Gene:", g, fill = TRUE)

             f <- dm_fitOneGeneManyGroups(y = counts[[g]], ngroups = ngroups, lgroups = lgroups, igroups = igroups, gamma0 = gamma0[g], prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)
             
             
             return(f)
             
           }, BPPARAM = BPPARAM))
           
           names(ff) <- gene_list
           
           pi <- MatrixList(lapply(ff, function(f) f$pi))
           stats <- do.call(rbind, lapply(ff, function(f) f$stats))
           
           cat("Took ", time["elapsed"], " seconds.\n")
           
           return(new("dmFit", proportions = pi, statistics = DataFrame(stats, row.names = rownames(stats))))
           
         },
         
         null = {
           
           cat("Fitting null model.. \n")
           group <- factor(rep("null", length(samples$group)))
           ngroups <- 1
           lgroups <- "null"
           igroups <- list(null = 1:length(samples$group))
           
           
           time <- system.time(ff <- BiocParallel::bplapply(gene_list, function(g){  
             # g = "ENSG00000135778"
             # cat("Gene:", g, fill = TRUE)
             
             f <- dm_fitOneGeneManyGroups(y = counts[[g]], ngroups = ngroups, lgroups = lgroups, igroups = igroups, gamma0 = gamma0[g], prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)
             
             return(f)
             
           }, BPPARAM = BPPARAM))
           
           names(ff) <- gene_list
           
           pi <- MatrixList(lapply(ff, function(f) f$pi))
           stats <- do.call(rbind, lapply(ff, function(f) f$stats))
           
           cat("Took ", time["elapsed"], " seconds.\n")
           
           return(new("dmFit", proportions = pi, statistics = DataFrame(stats, row.names = rownames(stats))))
           
           
         })
  
}


