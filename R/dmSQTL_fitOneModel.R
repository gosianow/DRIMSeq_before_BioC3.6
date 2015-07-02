##############################################################################
# multiple group fitting 
##############################################################################

dmSQTL_fitOneModel <- function(data, dispersion, model = c("full", "null")[1], prop_mode=c("constrOptim", "constrOptimG", "FisherScoring")[2], prop_tol = 1e-12, verbose=FALSE, BPPARAM = MulticoreParam(workers=1)){
  
  gene_list <- names(data$counts)
  
  if(length(dispersion == 1)){
    gamma0 <- lapply(data$genotypes, function(g){
      dispersion <- rep(dispersion, nrow(g))
      names(dispersion) <- rownames(g)
      return(dispersion)
    })
  } else {
    gamma0 <- dispersion
  }
  
  
  switch(model, 
         
         full={
           
           cat("Fitting full model.. \n")
           
           time <- system.time(fit <- bplapply(gene_list, function(g){
             # g = "ENSG00000163348.3"; y = data$counts[[g]]; snps = data$genotypes[[g]]
             
             y = data$counts[[g]]
             snps = data$genotypes[[g]]
             # snps <- snps[1:min(nrow(snps), 5), , drop = FALSE]
             
             f <- lapply(rownames(snps), function(i){
               # i = rownames(snps)[3]
               
               if(is.na(gamma0[[g]][i])) return(NULL)
               
               NAs <- is.na(snps[i, ]) | is.na(y[1, ])            
               yg <- y[, !NAs]             
               group <- snps[i, !NAs]
               group <- factor(group)
               ngroups <- nlevels(group)
               lgroups <- levels(group)
               nlibs <- length(group)
               
               igroups <- lapply(lgroups, function(gr){which(group == gr)})
               names(igroups) <- lgroups
               
               ff <- dm_fitOneGeneManyGroups(y = yg, ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
                                             gamma0 = gamma0[[g]][i], prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)  
               
               return(ff)
               
             })
             
             names(f) <- rownames(snps)
             
             return(f)
             
           }, BPPARAM = BPPARAM))
           
           cat("Took ", time["elapsed"], " seconds.\n")
           names(fit) <- gene_list
           return(fit)
           
         }, 
         
         null={
           
           cat("Fitting null model.. \n")
           
           time <- system.time(fit <- bplapply(gene_list, function(g){
             # g = gene_list[1]; y = data$counts[[g]]; snps = data$genotypes[[g]]
             
             y = data$counts[[g]]
             snps = data$genotypes[[g]]
             # snps <- snps[1:min(nrow(snps), 5), , drop = FALSE]
             
             groupg <- factor(rep("null", ncol(y)))
             ngroups <- 1
             lgroups <- "null"	 
             
             f <- lapply(rownames(snps), function(i){
               # i = rownames(snps)[3]
               
               if(is.na(gamma0[[g]][i]))
                 return(NULL)
               
               
               NAs <- is.na(snps[i, ]) | is.na(y[1, ])            
               yg <- y[, !NAs]             
               group <- groupg[!NAs]
               nlibs <- sum(!NAs)              
               igroups <- list(null = 1:nlibs)
               
               ff <- dm_fitOneGeneManyGroups(y = yg, ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
                                         gamma0 = gamma0[[g]][i], prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)
               
               return(ff)
               
             })
             
             names(f) <- rownames(snps)
            
             return(f)
             
           }, BPPARAM = BPPARAM))
           
           cat("Took ", time["elapsed"], " seconds.\n")
           names(fit) <- gene_list
           return(fit)
           
         })

}

