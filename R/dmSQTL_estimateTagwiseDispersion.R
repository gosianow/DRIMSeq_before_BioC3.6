##############################################################################
# calculate tagwise dispersions 
##############################################################################
# counts = x@counts; genotypes = x@genotypes; disp_adjust = TRUE; disp_mode = c("optimize", "optim", "constrOptim", "grid")[4]; disp_interval = c(0, 1e+5); disp_tol = 1e-08; disp_init = 100; disp_init_weirMoM = TRUE; disp_grid_length = 11; disp_grid_range = c(-10, 10); disp_moderation = c("none", "common", "trended")[1]; disp_prior_df = 0.1; disp_span = 0.2; prop_mode = c( "constrOptim", "constrOptimG", "FisherScoring")[2]; prop_tol = 1e-12; verbose = FALSE; BPPARAM = BiocParallel::MulticoreParam(workers = 10)

#' @importFrom stats optimize optim constrOptim complete.cases

dmSQTL_estimateTagwiseDispersion <- function(counts, genotypes, mean_expression, disp_adjust = TRUE, disp_mode = c("optimize", "optim", "constrOptim", "grid")[4], disp_interval = c(0, 1e+5), disp_tol = 1e-08,  disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 11, disp_grid_range = c(-10, 10), disp_moderation = c("none", "common", "trended")[1], disp_prior_df = 0.1, disp_span = 0.2, prop_mode = c("constrOptim", "constrOptimG", "FisherScoring")[2], prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
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
        
        not_nas <- complete.cases(loglik)        
        
        loglik <- loglik[not_nas, , drop = FALSE]
        
        
 
        if(disp_moderation != "none"){
          
          mean_expression <- rep(mean_expression, width(genotypes))[not_nas]
          
          ### Check where the grid is maximized 
          grid_max <- apply(loglik, 1, which.max)
          
          ### In the calculation of moderation, do not take into account genes that have dispersion on the top and bottom boundry of the grid (skipp 4 last grid points and 1 first grid point)
          not_boundry <- grid_max < (disp_grid_length - 3) & grid_max > 1
          boundry_last <- grid_max == disp_grid_length
          
          ### Calculate the span of the boundry loglikelihoods
          if(sum(boundry_last) > 1){
            loglik_span_boundry <- apply(loglik[boundry_last, , drop = FALSE], 1, function(x){max(x) - min(x)})
          }
          
          
          switch(disp_moderation, 
            
            common = {
              
              if(sum(not_boundry) == length(not_boundry)){
                moderation <- colMeans(loglik)
              }else{
                moderation <- colMeans(loglik[not_boundry, , drop = FALSE])
              }
              
              ### Estimate priorN
              ### Calculate the ratio between moderation lik span and lik span of boundry genes
              if(sum(boundry_last) > 1){

                moderation_span <- max(moderation) - min(moderation)
                
                span_ratio <- moderation_span / loglik_span_boundry

                # if(length(loglik_span_boundry) > 100){
                #   ### Do loess fitting if there is enough points
                #   df_priorN_loglog <- data.frame(priorN = log10(1/span_ratio), mean_expression = log10(mean_expression[boundry_last]))

                #   priorN_loess_loglog <- loess(priorN ~ mean_expression, df_priorN_loglog, control = loess.control(surface = "direct"))
                #   priorN_predict_loglog <- predict(priorN_loess_loglog, data.frame(mean_expression = log10(mean_expression)), se = FALSE)
                  
                #   priorN <- 10 ^ priorN_predict_loglog

                # }else{
                #   ### Otherwise, use median
                #   priorN <- quantile(1/span_ratio, 0.5)
                  
                # }
                
                priorN <- quantile(1/span_ratio, 0.5)

              }else{
                priorN <- disp_prior_df
              }
              
              if(length(priorN) == 1){
                message(paste0("! Using ", round(priorN, 4), " as a shrinkage factor !"))
                loglik <- sweep(loglik, 2, priorN * moderation, FUN = "+")
                }else{
                  message(paste0("! Using loess fit as a shrinkage factor !"))
                  loglik <- loglik + priorN * matrix(moderation, nrow = length(priorN), ncol = length(moderation), byrow = TRUE)
                }
              
              
            },
            
            trended = {
              
              
              
              if(sum(not_boundry) == length(not_boundry)){
                
                o <- order(mean_expression)
                oo <- order(o)
                width <- floor(disp_span * nrow(loglik))
                
                moderation <- edgeR::movingAverageByCol(loglik[o,], width = width)[oo,]
                
              }else{
                
                ### Use non boundry genes for calculating the moderation
                mean_expression_not_boundry <- mean_expression[not_boundry]
                loglik_not_boundry <- loglik[not_boundry, , drop = FALSE]
                
                o <- order(mean_expression_not_boundry)
                oo <- order(o)
                
                width <- floor(disp_span * nrow(loglik_not_boundry))
                
                moderation_not_boundry <- edgeR::movingAverageByCol(loglik_not_boundry[o, , drop = FALSE], width = width)[oo, , drop = FALSE]
                
                ### Fill in moderation values for the boundy genes
                moderation <- matrix(NA, nrow = nrow(loglik), ncol = ncol(loglik))
                
                moderation[not_boundry, ] <- moderation_not_boundry
                
                o <- order(mean_expression)
                oo <- order(o)
                
                moderation <- moderation[o, , drop = FALSE]
                not_boundry <- not_boundry[o]
                
                ### Last value in not_boundry must be TRUE
                if(not_boundry[length(not_boundry)] == FALSE){
                  
                  last_true <- max(which(not_boundry))
                  moderation[length(not_boundry), ] <- moderation[last_true, ]
                  
                  not_boundry[length(not_boundry)] <- TRUE
                  
                }
                
                not_boundry_diff <- diff(not_boundry, lag = 1)
                
                not_boundry_cumsum <- cumsum(not_boundry)
                
                ### Values used for filling in the boundry NAs - swith from FALSE to TRUE
                replacement_indx <- which(not_boundry_diff == 1) + 1
                
                replaced_indx <- which(!not_boundry)
                
                replaced_freq <- as.numeric(table(not_boundry_cumsum[replaced_indx]))
                
                moderation_boundry  <- moderation[rep(replacement_indx, times = replaced_freq), , drop = FALSE]
                
                moderation[!not_boundry, ] <- moderation_boundry
                
                moderation <- moderation[oo, , drop = FALSE]
                
              }
              
              ### Estimate priorN
              ### Calculate the ratio between moderation lik span and lik span of boundry genes
              if(sum(boundry_last) > 1){
                
                moderation_span_boundry <- apply(moderation[boundry_last, , drop = FALSE], 1, function(x){max(x) - min(x)})
                
                span_ratio <- moderation_span_boundry / loglik_span_boundry
                
                # if(length(loglik_span_boundry) > 100){
                #   ### Do loess fitting if there is enough points
                #   df_priorN_loglog <- data.frame(priorN = log10(1/span_ratio), mean_expression = log10(mean_expression[boundry_last]))

                #   priorN_loess_loglog <- loess(priorN ~ mean_expression, df_priorN_loglog, control = loess.control(surface = "direct"))
                #   priorN_predict_loglog <- predict(priorN_loess_loglog, data.frame(mean_expression = log10(mean_expression)), se = FALSE)
                  
                #   priorN <- 10 ^ priorN_predict_loglog

                # }else{
                #   ### Otherwise, use median
                #   priorN <- quantile(1/span_ratio, 0.5)
                  
                # }
                
                priorN <- quantile(1/span_ratio, 0.5)
                
              }else{
                priorN <- disp_prior_df
              }
              
              if(length(priorN) == 1){
                message(paste0("! Using ", round(priorN, 4), " as a shrinkage factor !"))
                loglik <- loglik + priorN * moderation
                }else{
                  message(paste0("! Using loess fit as a shrinkage factor !"))
                  loglik <- loglik + priorN * moderation
                }
              
            }
          )
          
        }
        
        out <- edgeR::maximizeInterpolant(splinePts, loglik)
        
        #### set NA for genes that tagwise disp could not be calculated            
        dispersion <- rep(NA, length(not_nas))
        names(dispersion) <- rownames(genotypes@unlistData)
        dispersion[not_nas] <- disp_init * 2^out
        dispersion <- relist(dispersion, genotypes@partitioning)
        
      }))
  
  if(verbose) cat("Took ", time["elapsed"], " seconds.\n")
  if(verbose) cat("*** Genewise dispersion: ", head(unlist(dispersion[1:6])), "... \n")
  
  return(dispersion)
}






