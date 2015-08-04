##############################################################################
# calculate tagwise dispersions 
##############################################################################

dmDS_estimateTagwiseDispersion <- function(counts, samples, mean_expression, disp_adjust = TRUE, disp_mode = c("optimize", "optim", "constrOptim", "grid")[4], disp_interval = c(0, 1e+5), disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = c("none", "common", "trended")[1], disp_prior_df = 10, disp_span = 0.3, prop_mode = c( "constrOptim", "constrOptimG", "FisherScoring")[2], prop_tol = 1e-12, verbose = FALSE, BPPARAM = MulticoreParam(workers=1)){
  
  gene_list <- names(counts)
  
  group <- samples$group
  ngroups <- nlevels(group)
  lgroups <- levels(group)
  nlibs <- length(group)
  
  igroups <- lapply(lgroups, function(gr){which(group == gr)})
  names(igroups) <- lgroups
  
  
  ### Find optimized dispersion
  cat("Estimating tagwise dispersion.. \n")
  time <- system.time(
    switch(
      disp_mode, 
      
      optimize={
        
        disp_list <- bplapply(gene_list, function(g){
          # g = 1
          #     print(g)
          
          ### return NA if gene has 1 exon or observations in one sample in group (anyway this gene would not be fitted by dmFit)
          if(is.na(dm_profileLikTagwise(gamma0 = disp_interval[1] + (1-(sqrt(5) - 1)/2)*(disp_interval[2]-disp_interval[1]), y = counts[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)))
            return(NA) 
          
          optimum <- optimize(f = dm_profileLikTagwise, interval = disp_interval,
                              y = counts[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, 
                              disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose,
                              maximum = TRUE, tol = disp_tol) 
          
          return(optimum$maximum)  
          
        }, BPPARAM = BPPARAM )
        
        
        names(disp_list) <- gene_list  
        dispersion <- unlist(disp_list)
        
      },
      
      
      optim={
        
        disp_list <- bplapply(gene_list, function(g){
          # g = 12
          if(verbose)
            cat("gene", g, "\n")
          
          ### return NA if gene has 1 exon or observations in one sample in group (anyway this gene would not be fitted by dmFit)
          if(is.na(dm_profileLikTagwise(gamma0 = disp_interval[1] + (1-(sqrt(5) - 1)/2)*(disp_interval[2]-disp_interval[1]), y = counts[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)))
            return(NA) 
          
          
          if(disp_init_weirMoM)
            disp_init <- dm_weirMoM(y = counts[[g]], se=FALSE)
          #              print(disp_init)
          
          
          optimum <- NA
          
          try( optimum <- optim(par = disp_init, fn = dm_profileLikTagwise, gr = NULL, 
                                y = counts[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, 
                                disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose,
                                method = "L-BFGS-B", lower = 1e-2, upper = 1e+10, control = list(fnscale = -1, factr = disp_tol))$par , silent = TRUE)
          
          
          
          return(optimum)  
          
        }, BPPARAM = BPPARAM)
        
        
        names(disp_list) <- gene_list  
        dispersion <- unlist(disp_list)
        
      },
      
      
      constrOptim={
        
        disp_list <- bplapply(gene_list, function(g){
          # g = 1
          #     print(g)
          
          ### return NA if gene has 1 exon or observations in one sample in group (anyway this gene would not be fitted by dmFit)
          if(is.null(dm_profileLikTagwise(gamma0 = disp_interval[1] + (1-(sqrt(5) - 1)/2)*(disp_interval[2]-disp_interval[1]), y = counts[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)))
            return(NA) 
          
          ui <- 1
          ci <- 1e-8
          
          if(disp_init_weirMoM)
            disp_init <- dm_weirMoM(y = counts[[g]], se=FALSE)
          
          
          optimum <- constrOptim(theta = disp_init, dm_profileLikTagwise, grad = NULL, method = "Nelder-Mead",
                                 ui=ui, ci=ci, control=list(fnscale = -1, reltol = disp_tol), 
                                 y = counts[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)
          
          
          return(optimum$par) 
          
        }, BPPARAM = BPPARAM )
        
        names(disp_list) <- gene_list  
        dispersion <- unlist(disp_list)
        
      },
      
      
      grid={
        
        ### genrate spline dispersion
        splinePts <- seq(from = disp_grid_range[1], to = disp_grid_range[2], length = disp_grid_length)
        splineDisp <- disp_init * 2^splinePts
        
        ### calculate the likelihood for each gene at the spline dispersion points
        
        loglikL <- bplapply(gene_list, function(g){
          # g = 1
          
          ll <- numeric(disp_grid_length)
          
          for(i in seq(disp_grid_length)){
            # i = 10 
            
            out <- dm_profileLikTagwise(gamma0 = splineDisp[i], y = counts[[g]], ngroups = ngroups, lgroups = lgroups, igroups = igroups, disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)
            
            if(is.na(out)){
              ll <- rep(NA, disp_grid_length)
              break
            }
            
            ll[i] <- out
            
          }
          
          return(ll)
          
        }, BPPARAM = BPPARAM)
        
        loglik <- do.call(rbind, loglikL)
        
        NAs <- complete.cases(loglik)        
        
        loglik <- loglik[NAs, , drop = FALSE]
        
        
        if(nrow(loglik) == 0){
          dispersion <- rep(NA, length(gene_list))
          names(dispersion) <- gene_list
          cat("** Tagwise dispersion: ", head(dispersion), "... \n")
          return(dispersion)
        }
        
        
        if(disp_moderation != "none"){
          
          nlibs <- length(group)
          
          priorN <- disp_prior_df/(nlibs - ngroups) ### analogy to edataR
          
          
          switch(
            disp_moderation, 
            
            common={
              
              moderation <- colMeans(loglik)
              
              loglik <- sweep(loglik, 1, priorN * moderation, FUN = "+")
              
            },
            
            trended={
              
              o <- order(mean_expression[NAs])
              oo <- order(o)
              width <- floor(disp_span * nrow(loglik))
              
              moderation <- edgeR::movingAverageByCol(loglik[o,], width=width)[oo,]
              
              loglik <- loglik + priorN * moderation ### like in edataR estimateTagwiseDisp
              # loglik <- (loglik + priorN * moderation)/(1 + priorN) ### like in edataR dispCoxReidInterpolateTagwise
              
            }
          )
          
        }
        
        
        out <- edgeR::maximizeInterpolant(splinePts, loglik)
        
        
        #### set NA for genes that tagwise disp could not be calculated 
        dispersion <- rep(NA, length(gene_list))
        names(dispersion) <- gene_list
        dispersion[NAs] <- disp_init * 2^out
        
        
      }))
  
  cat("Took ", time["elapsed"], " seconds.\n")
  cat("** Tagwise dispersion: ", head(dispersion), "... \n")
  
  return(dispersion)
  
}






