##############################################################################
# calculate tagwise dispersions 
##############################################################################

### Functions for bplapply

dmDS_grid_dm_profileLikTagwise <- function(g, seq_disp_grid_length, splineDisp, counts, ngroups, lgroups, igroups,  disp_adjust, prop_mode, prop_tol, verbose, disp_grid_length){
  # g = 1237
  
  ll <- numeric(disp_grid_length)
  
  for(i in seq_disp_grid_length){
    # i = 1
    
    out <- dm_profileLikTagwise(gamma0 = splineDisp[i], y = counts[[g]], 
      ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
      disp_adjust = disp_adjust, prop_mode = prop_mode, 
      prop_tol = prop_tol, verbose = verbose)
    
    if(is.na(out)){
      ll <- rep(NA, disp_grid_length)
      break
    }
    
    ll[i] <- out
    
  }
  
  return(ll)
  
}





# counts = x@counts; samples = x@samples; mean_expression = x@mean_expression; disp_adjust = TRUE; disp_mode = c("optimize", "optim", "constrOptim", "grid")[4]; disp_interval = c(0, 1e+5); disp_tol = 1e-08; disp_init = 100; disp_init_weirMoM = TRUE; disp_grid_length = 21; disp_grid_range = c(-10, 10); disp_moderation = c("none", "common", "trended")[1]; disp_prior_df = 0.1; disp_span = 0.2; prop_mode = c( "constrOptim", "constrOptimG", "FisherScoring")[2]; prop_tol = 1e-12; verbose = FALSE; BPPARAM = BiocParallel::MulticoreParam(workers = 10)

#' @importFrom stats optimize optim constrOptim complete.cases

dmDS_estimateTagwiseDispersion <- function(counts, samples, mean_expression, disp_adjust = TRUE, disp_mode = c("optimize", "optim", "constrOptim", "grid")[4], disp_interval = c(0, 1e+5), disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = c("none", "common", "trended")[1], disp_prior_df = 0, disp_span = 0.1, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  inds <- 1:length(counts)
  
  group <- samples$group
  ngroups <- nlevels(group)
  lgroups <- levels(group)
  nlibs <- length(group)
  
  igroups <- lapply(lgroups, function(gr){which(group == gr)})
  names(igroups) <- lgroups
  
  
  ### Find optimized dispersion
  if(verbose) cat("* Estimating genewise dispersion.. \n")
  time <- system.time(
    switch(
      disp_mode, 
      
      optimize={
        
        disp_list <- BiocParallel::bplapply(inds, function(g){
          # g = 1
          
          ### return NA if gene has 1 exon or observations in one sample in group (anyway this gene would not be fitted by dmFit)
          if(is.na(dm_profileLikTagwise(gamma0 = disp_interval[1] + (1-(sqrt(5) - 1)/2)*(disp_interval[2]-disp_interval[1]), y = counts[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)))
            return(NA) 
          
          optimum <- optimize(f = dm_profileLikTagwise, interval = disp_interval,
            y = counts[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, 
            disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose,
            maximum = TRUE, tol = disp_tol) 
          
          return(optimum$maximum)  
          
        }, BPPARAM = BPPARAM )
        
        names(disp_list) <- names(counts)  
        dispersion <- unlist(disp_list)
        
      },
      
      
      optim={
        
        disp_list <- BiocParallel::bplapply(inds, function(g){
          # g = 12
          
          ### return NA if gene has 1 exon or observations in one sample in group (anyway this gene would not be fitted by dmFit)
          if(is.na(dm_profileLikTagwise(gamma0 = disp_interval[1] + (1-(sqrt(5) - 1)/2)*(disp_interval[2]-disp_interval[1]), y = counts[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)))
            return(NA) 
          
          
          if(disp_init_weirMoM){
            disp_init_tmp <- dm_weirMoM(y = counts[[g]], se=FALSE)
            if(is.na(disp_init_tmp))
              disp_init_tmp <- disp_init
          }else{
            disp_init_tmp <- disp_init
          }
          
          
          optimum <- NA
          
          try( optimum <- optim(par = disp_init_tmp, fn = dm_profileLikTagwise, gr = NULL, 
            y = counts[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, 
            disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose,
            method = "L-BFGS-B", lower = 1e-2, upper = 1e+10, control = list(fnscale = -1, factr = disp_tol))$par , silent = TRUE)
          
          
          
          return(optimum)  
          
        }, BPPARAM = BPPARAM)
        
        
        names(disp_list) <- names(counts)  
        dispersion <- unlist(disp_list)
        
      },
      
      
      constrOptim={
        
        disp_list <- BiocParallel::bplapply(inds, function(g){
          # g = 1
          
          ### return NA if gene has 1 exon or observations in one sample in group (anyway this gene would not be fitted by dmFit)
          if(is.na(dm_profileLikTagwise(gamma0 = disp_interval[1] + (1-(sqrt(5) - 1)/2)*(disp_interval[2]-disp_interval[1]), y = counts[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)))
            return(NA) 
          
          ui <- 1
          ci <- 1e-8
          
          if(disp_init_weirMoM){
            disp_init_tmp <- dm_weirMoM(y = counts[[g]], se=FALSE)
            if(is.na(disp_init_tmp))
              disp_init_tmp <- disp_init
          }else{
            disp_init_tmp <- disp_init
          }
          
          
          optimum <- constrOptim(theta = disp_init_tmp, dm_profileLikTagwise, grad = NULL, method = "Nelder-Mead",
            ui=ui, ci=ci, control=list(fnscale = -1, reltol = disp_tol), 
            y = counts[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)
          
          
          return(optimum$par) 
          
        }, BPPARAM = BPPARAM )
        
        names(disp_list) <- names(counts)  
        dispersion <- unlist(disp_list)
        
      },
      
      
      grid={
        
        ### Standard grid from edgeR
        splinePts <- seq(from = disp_grid_range[1], to = disp_grid_range[2], length = disp_grid_length)
        
        ### More dense grid toward the common dispersion 
        # splinePts_uni <- sort(unique(c(0, seq(from = disp_grid_range[1], to = disp_grid_range[2], length = disp_grid_length))))
        
        # nr_positive_splitting <- sum(sign(splinePts_uni) == 1)
        # nr_negative_splitting <- sum(sign(splinePts_uni) == -1)
        
        # max_splitting <- max(nr_positive_splitting, nr_negative_splitting)
        # min_splitting <- min(nr_positive_splitting, nr_negative_splitting)
        
        # if(nr_positive_splitting == max_splitting)
        #   nr_splitting <- c((max_splitting - min_splitting + 1):max_splitting, max_splitting:1) + 2
        
        # if(nr_negative_splitting == max_splitting)
        #   nr_splitting <- c(1:max_splitting, max_splitting:(max_splitting - min_splitting + 1)) + 2
        
        # splinePts <- lapply(1:(length(splinePts_uni) - 1), function(i){
        
        #   seq(from = splinePts_uni[i], to = splinePts_uni[i + 1], length = nr_splitting[i])
        
        # })
        
        # splinePts <- sort(unique(unlist(splinePts)))
        
        
        
        disp_grid_length <- length(splinePts)
        splineDisp <- disp_init * 2^splinePts
        
        
        ### Calculate the likelihood for each gene at the spline dispersion points
        seq_disp_grid_length <- seq(disp_grid_length)
        
        loglikL <- BiocParallel::bplapply(inds, dmDS_grid_dm_profileLikTagwise, 
          seq_disp_grid_length = seq_disp_grid_length, splineDisp = splineDisp,
          counts = counts, 
          ngroups = ngroups, lgroups = lgroups, igroups = igroups,
          disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, 
          verbose = verbose, disp_grid_length = disp_grid_length, 
          BPPARAM = BPPARAM)
        
        
        loglik <- do.call(rbind, loglikL)
        not_nas <- complete.cases(loglik)        
        loglik <- loglik[not_nas, , drop = FALSE]
        
        if(nrow(loglik) == 0){
          dispersion <- rep(NA, length(inds))
          names(dispersion) <- names(counts)
          if(verbose) cat("*** Genewise dispersion: ", head(dispersion), "... \n")
          return(dispersion)
        }
        
        
        if(disp_moderation != "none"){
          
          mean_expression <- mean_expression[not_nas]
          
          loglik <- dm_profileLikModeration(loglik = loglik, mean_expression, disp_moderation = disp_moderation, disp_prior_df, disp_span)
          
          
        }
        
        
        out <- edgeR::maximizeInterpolant(splinePts, loglik)
        
        
        #### set NA for genes that tagwise disp could not be calculated 
        dispersion <- rep(NA, length(inds))
        names(dispersion) <- names(counts)
        dispersion[not_nas] <- disp_init * 2^out
        
        
      }))
  
  
  if(verbose) cat("Took ", time["elapsed"], " seconds.\n")
  if(verbose) cat("*** Genewise dispersion: ", head(dispersion), "... \n")
  
  return(dispersion)
  
}













