##############################################################################
# multiple group ffting 
##############################################################################
# counts = x@counts; samples = x@samples; dispersion = 38196.6; model = c("full", "null")[1]; prop_mode=c("constrOptim2", "constrOptim2G", "FisherScoring")[2]; prop_tol = 1e-12; verbose = FALSE; BPPARAM = BiocParallel::MulticoreParam(workers = 10)

dmDS_fitOneModel <- function(counts, samples, dispersion, model = c("full", "null")[1], prop_mode=c("constrOptim2", "constrOptim2G", "FisherScoring")[2], prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  inds <-  1:length(counts)
  
  if(length(dispersion) == 1){
    gamma0 <- rep(dispersion, length(inds))
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
           
           
           time <- system.time(ff <- BiocParallel::bplapply(inds, function(g){  
             # g = 1
             # cat("Gene:", g, fill = TRUE)

             f <- dm_fitOneGeneManyGroups(y = counts[[g]], ngroups = ngroups, lgroups = lgroups, igroups = igroups, gamma0 = gamma0[g], prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)
             
             
             return(f)
             
           }, BPPARAM = BPPARAM))
           
           names(ff) <- names(counts)
           
           stats <- do.call(rbind, lapply(ff, function(f) f[[2]])) ### stats
           rownames(stats) <- names(counts)
           
           fff <- MatrixList(lapply(ff, function(f) f[[1]]), metadata = stats) ### pi 
           
           cat("Took ", time["elapsed"], " seconds.\n")
           
           return(fff)
           
         },
         
         null = {
           
           cat("Fitting null model.. \n")
           group <- factor(rep("null", length(samples$group)))
           ngroups <- 1
           lgroups <- "null"
           igroups <- list(null = 1:length(samples$group))
           
           
           time <- system.time(ff <- BiocParallel::bplapply(inds, function(g){  
             # g = 1
             # cat("Gene:", g, fill = TRUE)
             
             f <- dm_fitOneGeneManyGroups(y = counts[[g]], ngroups = ngroups, lgroups = lgroups, igroups = igroups, gamma0 = gamma0[g], prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)
             
             return(f)
             
           }, BPPARAM = BPPARAM))
           
           names(ff) <- names(counts)
           
           stats <- do.call(rbind, lapply(ff, function(f) f[[2]])) ### stats
           rownames(stats) <- names(counts)
           
           fff <- MatrixList(lapply(ff, function(f) f[[1]]), metadata = stats) ### pi 
           
           cat("Took ", time["elapsed"], " seconds.\n")
           
           return(fff)
           
           
         })
  
}


