##############################################################################
# calculate tagwise dispersions 
##############################################################################

# dge = dgeDM; adjustDisp = TRUE; modeDisp = modeDisp; intervalDisp = c(0, 1e+5); tolDisp = 1e-08;  initDisp = 10; initWeirMoMDisp = TRUE; gridLengthDisp = 15; gridRangeDisp = c(-7, 7); trendDisp = trendDisp; priorDfDisp = 10; spanDisp = 0.3; modeProp = modeProp; tolProp = 1e-12; verbose = FALSE; plot = FALSE; BPPARAM = BPPARAM

dmEstimateTagwiseDisp <- function(dge, adjustDisp = TRUE, modeDisp = c("optimize", "optim", "constrOptim", "grid")[2], intervalDisp = c(0, 1e+5), tolDisp = 1e-08,  initDisp = 10, initWeirMoMDisp = TRUE, gridLengthDisp = 15, gridRangeDisp = c(-7, 7), trendDisp = c("none", "commonDispersion", "trendedDispersion")[1], priorDfDisp = 10, spanDisp = 0.3, modeProp = c( "constrOptim2", "constrOptim2G", "FisherScoring")[2], tolProp = 1e-12, verbose = FALSE, plot = FALSE, BPPARAM = MulticoreParam(workers=1)){
  
  ### calculate mean expression of genes 
	cat("Calculating mean gene expression.. \n")
  time <- system.time(meanExpr <- unlist(bplapply(dge$counts, function(g){ mean(colSums(g), na.rm = TRUE) }, BPPARAM = BPPARAM)))
	cat("Took ", time["elapsed"], " seconds.\n")
  dge$meanExpr <- meanExpr
  geneList <- names(dge$counts)
  

  group <- factor(dge$samples$group)
  ngroups <- nlevels(group)
  lgroups <- levels(group)
  nlibs <- length(group)
  
  igroups <- list()
  for(gr in 1:ngroups){
    # gr=2
    igroups[[lgroups[gr]]] <- which(group == lgroups[gr])
    
  }
 
  
  ### Find optimized dispersion
	cat("Estimating tagwise dispersion.. \n")
  time <- system.time(switch(modeDisp, 
         
         optimize={
           
           dispList <- bplapply(geneList, function(g){
             # g = 1
             #     print(g)
             
             ### return NA if gene has 1 exon or observations in one sample in group (anyway this gene would not be fitted by dmFit)
             if(is.null(dmAdjustedProfileLikTG(gamma0 = intervalDisp[1] + (1-(sqrt(5) - 1)/2)*(intervalDisp[2]-intervalDisp[1]), y = dge$counts[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjustDisp = adjustDisp, modeProp = modeProp, tolProp = tolProp, verbose = verbose)))
							 		return(NA) 
             
             optimum <- optimize(f = dmAdjustedProfileLikTG, interval = intervalDisp,
                                 y = dge$counts[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, 
                                 adjustDisp = adjustDisp, modeProp = modeProp, tolProp = tolProp, verbose = verbose,
                                 maximum = TRUE, tol = tolDisp) 
             
             return(optimum$maximum)  
             
           }, BPPARAM = BPPARAM )
           
           
           names(dispList) <- geneList  
           dge$tagwiseDispersion <- unlist(dispList)
           
         },
         
         
         optim={
           
      dispList <- bplapply(geneList, function(g){
             # g = 12
						 if(verbose)
							 cat("gene", g, "\n")
						 
             ### return NA if gene has 1 exon or observations in one sample in group (anyway this gene would not be fitted by dmFit)
             if(is.null(dmAdjustedProfileLikTG(gamma0 = intervalDisp[1] + (1-(sqrt(5) - 1)/2)*(intervalDisp[2]-intervalDisp[1]), y = dge$counts[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjustDisp = adjustDisp, modeProp = modeProp, tolProp = tolProp, verbose = verbose)))
               return(NA) 
             
             
             if(initWeirMoMDisp)
               initDisp <- weirMoM(data = dge$counts[[g]], se=FALSE)
             #              print(initDisp)
             
             
             
             try( optimum <- optim(par = initDisp, fn = dmAdjustedProfileLikTG, gr = NULL, 
                          y = dge$counts[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, 
                                 adjustDisp = adjustDisp, modeProp = modeProp, tolProp = tolProp, verbose = verbose,
                          method = "L-BFGS-B", lower = 1e-2, upper = 1e+10, control = list(fnscale = -1, factr = tolDisp)) , silent = TRUE)
             
             #              print(out$par)
             
             
             return(optimum$par)  
             
           }, BPPARAM = BPPARAM)
           
           
           names(dispList) <- geneList  
           dge$tagwiseDispersion <- unlist(dispList)
                  
         },



         constrOptim={
                    
           dispList <- bplapply(geneList, function(g){
             # g = 1
             #     print(g)
             
             ### return NA if gene has 1 exon or observations in one sample in group (anyway this gene would not be fitted by dmFit)
             if(is.null(dmAdjustedProfileLikTG(gamma0 = intervalDisp[1] + (1-(sqrt(5) - 1)/2)*(intervalDisp[2]-intervalDisp[1]), y = dge$counts[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjustDisp = adjustDisp, modeProp = modeProp, tolProp = tolProp, verbose = verbose)))
               return(NA) 
             
             ui <- 1
             ci <- 1e-8
             
             if(initWeirMoMDisp)
               initDisp <- weirMoM(data = dge$counts[[g]], se=FALSE)
             
             
             optimum <- constrOptim(theta = initDisp, dmAdjustedProfileLikTG, grad = NULL, method = "Nelder-Mead",
                                ui=ui, ci=ci, control=list(fnscale = -1, reltol = tolDisp), 
                                y = dge$counts[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjustDisp = adjustDisp, modeProp = modeProp, tolProp = tolProp, verbose = verbose)
             
             
             return(optimum$par) 
             
           }, BPPARAM = BPPARAM )
           
           names(dispList) <- geneList  
           dge$tagwiseDispersion <- unlist(dispList)
           
         },
      
      
          grid={
        
            ### genrate spline dispersion
            splinePts <- seq(from = gridRangeDisp[1], to = gridRangeDisp[2], length = gridLengthDisp)
            splineDisp <- dge$commonDispersion * 2^splinePts
            
            ### calculate the likelihood for each gene at the spline dispersion points
            
            loglik0L <- bplapply(geneList, function(g){
              # g = 1
              
              ll <- numeric(gridLengthDisp)
              
              for(i in seq(gridLengthDisp)){
                # i = 10 
                out <- dmAdjustedProfileLikTG(gamma0 = splineDisp[i], y = dge$counts[[g]], ngroups = ngroups, lgroups = lgroups, igroups = igroups, adjustDisp = adjustDisp, modeProp = modeProp, tolProp = tolProp, verbose = verbose)
                if(is.null(out))
                  return(NULL)
                
                ll[i] <- out
                
              }
              
              return(ll)
              
            }, BPPARAM = BPPARAM)
            
            names(loglik0L) <- geneList            
            loglik0 <- do.call(rbind, loglik0L)
            
						
						if(is.null(loglik0)){
	            dge$tagwiseDispersion <- rep(NA, length(geneList))
	            names(dge$tagwiseDispersion) <- geneList
							cat("** Tagwise dispersion: ", head(dge$tagwiseDispersion), "... \n")
							return(dge)
						}
						
						
            genesComplete <- rownames(loglik0)
            loglik <- loglik0 
            
            if(trendDisp != "none"){
  
              switch(trendDisp, 
                     commonDispersion={

                       moderation <- matrix(colMeans(loglik0), nrow(loglik0), gridLengthDisp, byrow=TRUE)
                       
                     },
                     
                     trendedDispersion={
                       
                       o <- order(dge$meanExpr[genesComplete])
                       oo <- order(o)
                       width <- floor(spanDisp * nrow(loglik0))
                       
                       moderation <- edgeR::movingAverageByCol(loglik0[o,], width=width)[oo,]
                       
                     })
              
              rownames(moderation) <- genesComplete
              nlibs <- length(group)
							
              priorN <- priorDfDisp/(nlibs - ngroups) ### analogy to edgeR
              
              loglik <- loglik0 + priorN * moderation ### like in edgeR estimateTagwiseDisp
#               loglik <- (loglik0 + priorN * moderation)/(1 + priorN) ### like in edgeR dispCoxReidInterpolateTagwise

if(plot){
  
  dge$plotSplineDisp <- splineDisp
  dge$plotLoglik0 <- loglik0
  dge$plotModeration <- moderation
  dge$plotPriorN <- priorN
  dge$plotLoglik <- loglik
  
}
              
            }
              

            out <- edgeR::maximizeInterpolant(splinePts, loglik)
            
            names(out) <- genesComplete
            
if(plot){
  dge$plotOutDisp <- dge$commonDispersion * 2^out
}
#### set common dispersion for genes that tagwise disp could not be calculated 
            dge$tagwiseDispersion <- rep(NA, length(geneList))
            names(dge$tagwiseDispersion) <- geneList
            dge$tagwiseDispersion[genesComplete] <- dge$commonDispersion * 2^out
  
        
      }))
  
  cat("Took ", time["elapsed"], " seconds.\n")
  cat("** Tagwise dispersion: ", head(dge$tagwiseDispersion), "... \n")
  
  return(dge)
  
}






