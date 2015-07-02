##############################################################################
# calculate tagwise dispersions 
##############################################################################

dmDS_estimateTagwiseDispersion <- function(data, disp_adjust = TRUE, disp_mode = c("optimize", "optim", "constrOptim", "grid")[4], disp_interval = c(0, 1e+5), disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = c("none", "common", "trended")[1], disp_prior_df = 10, disp_span = 0.3, prop_mode = c( "constrOptim", "constrOptimG", "FisherScoring")[2], prop_tol = 1e-12, verbose = FALSE, BPPARAM = MulticoreParam(workers=1)){
  
  gene_list <- names(data$counts)
  
  group <- data$samples$group
  ngroups <- nlevels(group)
  lgroups <- levels(group)
  nlibs <- length(group)
  
  igroups <- lapply(lgroups, function(gr){which(group == gr)})
  names(igroups) <- lgroups
 
  
  ### Find optimized dispersion
	cat("Estimating tagwise dispersion.. \n")
  time <- system.time(switch(disp_mode, 
         
         optimize={
           
           disp_list <- bplapply(gene_list, function(g){
             # g = 1
             #     print(g)
             
             ### return NA if gene has 1 exon or observations in one sample in group (anyway this gene would not be fitted by dmFit)
             if(is.null(dm_profileLikTagwise(gamma0 = disp_interval[1] + (1-(sqrt(5) - 1)/2)*(disp_interval[2]-disp_interval[1]), y = data$counts[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)))
							 		return(NA) 
             
             optimum <- optimize(f = dm_profileLikTagwise, interval = disp_interval,
                                 y = data$counts[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, 
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
             if(is.null(dm_profileLikTagwise(gamma0 = disp_interval[1] + (1-(sqrt(5) - 1)/2)*(disp_interval[2]-disp_interval[1]), y = data$counts[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)))
               return(NA) 
             
             
             if(disp_init_weirMoM)
               disp_init <- dm_weirMoM(y = data$counts[[g]], se=FALSE)
             #              print(disp_init)
             
             
             
             try( optimum <- optim(par = disp_init, fn = dm_profileLikTagwise, gr = NULL, 
                          y = data$counts[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, 
                                 disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose,
                          method = "L-BFGS-B", lower = 1e-2, upper = 1e+10, control = list(fnscale = -1, factr = disp_tol)) , silent = TRUE)
             
             #              print(out$par)
             
             
             return(optimum$par)  
             
           }, BPPARAM = BPPARAM)
           
           
           names(disp_list) <- gene_list  
           dispersion <- unlist(disp_list)
                  
         },


         constrOptim={
                    
           disp_list <- bplapply(gene_list, function(g){
             # g = 1
             #     print(g)
             
             ### return NA if gene has 1 exon or observations in one sample in group (anyway this gene would not be fitted by dmFit)
             if(is.null(dm_profileLikTagwise(gamma0 = disp_interval[1] + (1-(sqrt(5) - 1)/2)*(disp_interval[2]-disp_interval[1]), y = data$counts[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)))
               return(NA) 
             
             ui <- 1
             ci <- 1e-8
             
             if(disp_init_weirMoM)
               disp_init <- dm_weirMoM(y = data$counts[[g]], se=FALSE)
             
             
             optimum <- constrOptim(theta = disp_init, dm_profileLikTagwise, grad = NULL, method = "Nelder-Mead",
                                ui=ui, ci=ci, control=list(fnscale = -1, reltol = disp_tol), 
                                y = data$counts[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)
             
             
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
            
            loglik0L <- bplapply(gene_list, function(g){
              # g = 1
              
              ll <- numeric(disp_grid_length)
              
              for(i in seq(disp_grid_length)){
                # i = 10 
                out <- dm_profileLikTagwise(gamma0 = splineDisp[i], y = data$counts[[g]], ngroups = ngroups, lgroups = lgroups, igroups = igroups, disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)
                if(is.null(out))
                  return(NULL)
                
                ll[i] <- out
                
              }
              
              return(ll)
              
            }, BPPARAM = BPPARAM)
            
            names(loglik0L) <- gene_list            
            loglik0 <- do.call(rbind, loglik0L)
            
						
						if(is.null(loglik0)){
	            dispersion <- rep(NA, length(gene_list))
	            names(dispersion) <- gene_list
							cat("** Tagwise dispersion: ", head(dispersion), "... \n")
							return(dispersion)
						}
						
						
            genesComplete <- rownames(loglik0)
            loglik <- loglik0 
            
            if(disp_moderation != "none"){
  
              switch(disp_moderation, 
                     common={

                       moderation <- matrix(colMeans(loglik0), nrow(loglik0), disp_grid_length, byrow=TRUE)
                       
                     },
                     
                     trended={
                       
                       o <- order(mean_expression[genesComplete])
                       oo <- order(o)
                       width <- floor(disp_span * nrow(loglik0))
                       
                       moderation <- edgeR::movingAverageByCol(loglik0[o,], width=width)[oo,]
                       
                     })
              
              rownames(moderation) <- genesComplete
              nlibs <- length(group)
							
              priorN <- disp_prior_df/(nlibs - ngroups) ### analogy to edataR
              
              loglik <- loglik0 + priorN * moderation ### like in edataR estimateTagwiseDisp
#               loglik <- (loglik0 + priorN * moderation)/(1 + priorN) ### like in edataR dispCoxReidInterpolateTagwise

# if(plot){
#
#   data$plotSplineDisp <- splineDisp
#   data$plotLoglik0 <- loglik0
#   data$plotModeration <- moderation
#   data$plotPriorN <- priorN
#   data$plotLoglik <- loglik
#
# }
              
            }
              

            out <- edgeR::maximizeInterpolant(splinePts, loglik)
            
            names(out) <- genesComplete
            
# if(plot){
#   data$plotOutDisp <- data$commonDispersion * 2^out
# }
#### set NA for genes that tagwise disp could not be calculated 
            dispersion <- rep(NA, length(gene_list))
            names(dispersion) <- gene_list
            dispersion[genesComplete] <- disp_init * 2^out
  
        
      }))
  
  cat("Took ", time["elapsed"], " seconds.\n")
  cat("** Tagwise dispersion: ", head(dispersion), "... \n")
	
  return(dispersion)
  
}






