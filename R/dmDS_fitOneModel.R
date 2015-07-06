##############################################################################
# multiple group fitting 
##############################################################################


dmDS_fitOneModel <- function(data, dispersion, model = c("full", "null")[1], prop_mode=c("constrOptim2", "constrOptim2G", "FisherScoring")[2], prop_tol = 1e-12, verbose = FALSE, BPPARAM = MulticoreParam(workers=1)){
  
  gene_list <-  names(data@counts)
	
  if(length(dispersion) == 1){
    gamma0 <- rep(dispersion, length(gene_list))
		names(gamma0) <- gene_list
  } else {
    gamma0 <- dispersion
  }
  
  
  switch(model, 
         
         full = {
           cat("Fitting full model.. \n")
           group <- data@samples$group
           ngroups <- nlevels(group)
           lgroups <- levels(group)
           
           igroups <- lapply(lgroups, function(gr){which(group == gr)})
           names(igroups) <- lgroups
           
           
           time <- system.time(fit <- bplapply(gene_list, function(g){  
             # g = "ENSG00000135778"
             # cat("Gene:", g, fill = TRUE)
             
             if(is.na(gamma0[g]))
               return(NULL)
             
             f <- dm_fitOneGeneManyGroups(y = data@counts[[g]], ngroups = ngroups, lgroups = lgroups, igroups = igroups, gamma0 = gamma0[g], prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)
             
             return(f)
             
           }, BPPARAM = BPPARAM))
           
           names(fit) <- gene_list
           
           cat("Took ", time["elapsed"], " seconds.\n")
           
           return(fit)
           
           
         },
         
         null = {
           
           cat("Fitting null model.. \n")
           group <- factor(rep("null", length(data@samples$group)))
           ngroups <- 1
           lgroups <- "null"
           igroups <- list(null = 1:length(data@samples$group))
           
           
           time <- system.time(fit <- bplapply(gene_list, function(g){  
             # g = "ENSG00000135778"
             # cat("Gene:", g, fill = TRUE)
             
             if(is.na(gamma0[g]))
               return(NULL)
             
             f <- dm_fitOneGeneManyGroups(y = data@counts[[g]], ngroups = ngroups, lgroups = lgroups, igroups = igroups, gamma0 = gamma0[g], prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)
             
             return(f)
             
           }, BPPARAM = BPPARAM))
           
           names(fit) <- gene_list
           
           cat("Took ", time["elapsed"], " seconds.\n")
           
           return(fit)
           
         })
  
}


