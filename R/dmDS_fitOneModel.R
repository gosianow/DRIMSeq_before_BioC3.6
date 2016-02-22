##############################################################################
# multiple group ffting 
##############################################################################

dmDS_fitOneModel <- function(counts, samples, dispersion, model = "full",
  prop_mode = "constrOptim2G", prop_tol = 1e-12, verbose = FALSE, 
  BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  inds <-  1:length(counts)
  
  if(length(dispersion) == 1){
    gamma0 <- rep(dispersion, length(inds))
  } else {
    gamma0 <- dispersion
  }
  
  switch(model, 
    
    full = {
      
      if(verbose) message("* Fitting full model.. \n")
      
      group <- samples$group
      ngroups <- nlevels(group)
      lgroups <- levels(group)
      
      igroups <- lapply(lgroups, function(gr){which(group == gr)})
      names(igroups) <- lgroups
 
      time <- system.time(ff <- BiocParallel::bplapply(inds, function(g, counts, 
        ngroups, lgroups, igroups, gamma0, prop_mode, prop_tol, verbose){  
        # g = 1
        if(verbose >= 2)
          message(" Gene:", g)
        
        f <- dm_fitOneGeneManyGroups(y = counts[[g]], ngroups = ngroups, 
          lgroups = lgroups, igroups = igroups, gamma0 = gamma0[g], 
          prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)
        
        return(f)
        
      }, counts = counts, ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
        gamma0 = gamma0, prop_mode = prop_mode, prop_tol = prop_tol, 
        verbose = verbose, BPPARAM = BPPARAM))
      
      names(ff) <- names(counts)
      
      stats <- do.call(rbind, lapply(ff, function(f) f[[2]])) ### stats: liks
      rownames(stats) <- names(counts)
      
      fff <- MatrixList(lapply(ff, function(f) f[[1]]), metadata = stats) ### pis and liks
      
      if(verbose >= 2) message("\n")
      if(verbose) message("Took ", time["elapsed"], " seconds.\n")
      
      return(fff)
      
    },
    
    null = {
      
      if(verbose) message("* Fitting null model.. \n")
      
      time <- system.time(ff <- BiocParallel::bplapply(inds, function(g, counts, 
        gamma0, prop_mode, prop_tol, verbose){  
        # g = 1
        # message("Gene:", g)
        
        f <- dm_fitOneGeneOneGroup(y = counts[[g]], gamma0 = gamma0[g], 
          prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)
        
        f[[1]] <- matrix(f[[1]], dimnames = list(names(f[[1]]), "null")) ### convert pi into matrix
        
        return(f)
        
      }, counts = counts, gamma0 = gamma0, prop_mode = prop_mode, 
        prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM))
      
      names(ff) <- names(counts)
      
      stats <- do.call(rbind, lapply(ff, function(f) f[[2]])) ### stats: lik, df
      rownames(stats) <- names(counts)
      colnames(stats) <- c("lik", "df")
      
      fff <- MatrixList(lapply(ff, function(f) f[[1]]), metadata = stats) ### pi 
      
      if(verbose) message("Took ", time["elapsed"], " seconds.\n")
      
      return(fff)

    })
  
}


