##############################################################################
# calculate tagwise dispersions 
##############################################################################
# counts = x@counts; genotypes = x@genotypes; disp_adjust = TRUE; disp_mode = c("optimize", "optim", "constrOptim", "grid")[4]; disp_interval = c(0, 1e+5); disp_tol = 1e-08; disp_init = 100; disp_init_weirMoM = TRUE; disp_grid_length = 21; disp_grid_range = c(-10, 10); disp_moderation = c("none", "common", "trended")[1]; disp_prior_df = 10; disp_span = 0.3; prop_mode = c( "constrOptim", "constrOptimG", "FisherScoring")[2]; prop_tol = 1e-12; verbose = FALSE; BPPARAM = BiocParallel::MulticoreParam(workers = 10)

#' @importFrom stats optimize optim constrOptim complete.cases

dmSQTL_estimateTagwiseDispersion <- function(counts, genotypes, mean_expression, disp_adjust = TRUE, disp_mode = c("optimize", "optim", "constrOptim", "grid")[4], disp_interval = c(0, 1e+5), disp_tol = 1e-08,  disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = c("none", "common", "trended")[1], disp_prior_df = 10, disp_span = 0.3, prop_mode = c("constrOptim", "constrOptimG", "FisherScoring")[2], prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  inds <- 1:length(counts)
  if(verbose) cat("* Estimating genewise dispersion.. \n")
  
  time <- system.time( 
    switch(
      disp_mode, 
      
      optimize={
        
        disp_list <- BiocParallel::bplapply(inds, function(g){
          # g = 1
          y = counts[[g]]
          snps = genotypes[[g]]
          disp <- rep(NA, nrow(snps))
          names(disp) <- rownames(snps)
          
          for(i in 1:nrow(snps)){
            # i = 1
            
            NAs <- is.na(snps[i, ]) | is.na(y[1, ])            
            yg <- y[, !NAs]             
            group <- snps[i, !NAs]
            group <- factor(group)
            ngroups <- nlevels(group)
            lgroups <- levels(group)
            nlibs <- length(group)
            
            igroups <- lapply(lgroups, function(gr){which(group == gr)})
            names(igroups) <- lgroups
            
            gamma0 <- disp_interval[1] + (1-(sqrt(5) - 1)/2)*(disp_interval[2]-disp_interval[1])
            if(is.na(dm_profileLikTagwise(gamma0 = gamma0, y = yg, ngroups = ngroups, lgroups = lgroups, igroups=igroups, disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)))
              next
            
            optimum <- optimize(f = dm_profileLikTagwise, interval = disp_interval,
                                y = yg, ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
                                disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose,
                                maximum = TRUE, tol = disp_tol)
            
            disp[i] <- optimum$maximum       
            
          }
          return(disp)  
        }, BPPARAM = BPPARAM)
        
        names(disp_list) <- names(counts)
        dispersion <- disp_list
        
      },
      
      optim={
        
        disp_list <- BiocParallel::bplapply(inds, function(g){
          # g = 1
          y = counts[[g]]
          snps = genotypes[[g]]
          disp <- rep(NA, nrow(snps))
          names(disp) <- rownames(snps)
          
          for(i in 1:nrow(snps)){
            # i = 1
            
            NAs <- is.na(snps[i, ]) | is.na(y[1, ])            
            yg <- y[, !NAs]             
            group <- snps[i, !NAs]
            group <- factor(group)
            ngroups <- nlevels(group)
            lgroups <- levels(group)
            nlibs <- length(group)
            
            igroups <- lapply(lgroups, function(gr){which(group == gr)})
            names(igroups) <- lgroups
            
            gamma0 <- disp_init
            if(is.na(dm_profileLikTagwise(gamma0 = gamma0, y = yg, ngroups=ngroups, lgroups=lgroups, igroups=igroups, disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)))
              next
            
            if(disp_init_weirMoM){
              disp_init_tmp <- dm_weirMoM(y = counts[[g]], se=FALSE)
              if(is.na(disp_init_tmp))
                disp_init_tmp <- disp_init
            }else{
              disp_init_tmp <- disp_init
            }
            
            try( optimum <- optim(par = disp_init_tmp, fn = dm_profileLikTagwise, gr = NULL, 
                                  y = yg, ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
                                  disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose,
                                  method = "L-BFGS-B", lower = 1e-2, upper = 1e+10, control = list(fnscale = -1, factr = disp_tol)), silent = TRUE )
            
            disp[i] <- optimum$par       
            
          }
          return(disp)  
        }, BPPARAM = BPPARAM)
        
        names(disp_list) <- names(counts)
        dispersion <- disp_list
        
      },
      
      constrOptim={
        
        disp_list <- BiocParallel::bplapply(inds, function(g){
          # g = 1
          y = counts[[g]]
          snps = genotypes[[g]]
          disp <- rep(NA, nrow(snps))
          names(disp) <- rownames(snps)
          
          for(i in 1:nrow(snps)){
            # i = 1
            
            NAs <- is.na(snps[i, ]) | is.na(y[1, ])            
            yg <- y[, !NAs]             
            group <- snps[i, !NAs]
            group <- factor(group)
            ngroups <- nlevels(group)
            lgroups <- levels(group)
            nlibs <- length(group)
            
            igroups <- lapply(lgroups, function(gr){which(group == gr)})
            names(igroups) <- lgroups
            
            gamma0 <- disp_init
            if(is.na(dm_profileLikTagwise(gamma0 = gamma0, y = yg, ngroups=ngroups, lgroups=lgroups, igroups=igroups, disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)))
              next
            
            ui <- 1
            ci <- 1e-8
            
            if(disp_init_weirMoM){
              disp_init_tmp <- dm_weirMoM(y = counts[[g]], se=FALSE)
              if(is.na(disp_init_tmp))
                disp_init_tmp <- disp_init
            }else{
              disp_init_tmp <- disp_init
            }
            
            optimum <- constrOptim(theta = disp_init_tmp, dm_profileLikTagwise, grad = NULL, method = "Nelder-Mead", ui=ui, ci=ci, control = list(fnscale = -1, reltol = disp_tol),  y = yg, ngroups=ngroups, lgroups=lgroups, igroups=igroups, disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose )

            disp[i] <- optimum$par
            
          }
          return(disp)
        }, BPPARAM = BPPARAM)
        
        names(disp_list) <- names(counts)
        dispersion <- disp_list
        
      },
      
      grid={
        
        ### genrate spline dispersion
        splinePts <- seq(from = disp_grid_range[1], to = disp_grid_range[2], length = disp_grid_length)
        splineDisp <- disp_init * 2^splinePts
        
        ### calculate the likelihood for each gene at the spline dispersion points
        seq_disp_grid_length <- seq(disp_grid_length)
        
        loglikL <- BiocParallel::bplapply(inds, function(g){
          # g = 1
          
          y = counts[[g]]
          snps = genotypes[[g]]
          ll <- matrix(0, nrow(snps), disp_grid_length)
          
          
          for(i in 1:nrow(snps)){
            # i = 1
            
            NAs <- is.na(snps[i, ]) | is.na(y[1, ])            
            yg <- y[, !NAs]             
            group <- snps[i, !NAs]
            group <- factor(group)
            ngroups <- nlevels(group)
            lgroups <- levels(group)
            nlibs <- length(group)
            
            igroups <- lapply(lgroups, function(gr){which(group == gr)})
            names(igroups) <- lgroups
            
            for(j in seq_disp_grid_length){
              # j = 1 
              
              out <- dm_profileLikTagwise(gamma0 = splineDisp[j], y = yg, ngroups = ngroups, lgroups = lgroups, igroups = igroups, disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)
              
              if(is.na(out)){
                ll[i, ] <- NA
                break
              }
              
              ll[i, j] <- out
              
            } # j
          } # i
          
          return(ll)
          
        }, BPPARAM = BPPARAM)
        
        loglik <- do.call(rbind, loglikL)
        NAs <- complete.cases(loglik)        
        
        loglik <- loglik[NAs, , drop = FALSE]
        
        if(disp_moderation != "none"){
          
          # ### FIX IT!
          # nlibs <- ncol(snps)
          # ngroups <- 2
          # priorN <- disp_prior_df/(nlibs - ngroups) ### analogy to edgeR
          
          priorN <- disp_prior_df
          
          switch(
            disp_moderation, 
            
            common={
              
              moderation <- colMeans(loglik)
              
              loglik <- sweep(loglik, 2, priorN * moderation, FUN = "+")
              
            },
            
            trended={
              
              mean_expression <- rep(mean_expression, width(genotypes))[NAs]
              o <- order(mean_expression)
              oo <- order(o)
              width <- floor(disp_span * nrow(loglik))
              
              moderation <- edgeR::movingAverageByCol(loglik[o,], width = width)[oo,]
              
              loglik <- loglik + priorN * moderation ### like in edgeR estimateTagwiseDisp
              # loglik <- (loglik + priorN * moderation)/(1 + priorN) ### like in edgeR dispCoxReidInterpolateTagwise
              
            }
          )
          
        }
        
        out <- edgeR::maximizeInterpolant(splinePts, loglik)
        
        #### set NA for genes that tagwise disp could not be calculated            
        dispersion <- rep(NA, length(NAs))
        names(dispersion) <- rownames(genotypes@unlistData)
        dispersion[NAs] <- disp_init * 2^out
        dispersion <- relist(dispersion, genotypes@partitioning)
        
      }))
  
  if(verbose) cat("Took ", time["elapsed"], " seconds.\n")
  if(verbose) cat("*** Genewise dispersion: ", head(unlist(dispersion[1:6])), "... \n")
  
  return(dispersion)
}






