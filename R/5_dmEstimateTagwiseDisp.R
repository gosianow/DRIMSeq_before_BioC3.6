##############################################################################
# calculate tagwise dispersions 
##############################################################################

# group=NULL; adjust = TRUE; mode = "constrOptim2G"; epsilon = 1e-05; maxIte = 1000; modeDisp=c("optimize", "optim", "constrOptim", "grid")[4]; interval = c(0, 1e+5); tol = 1e-00;  initDisp = 10; initWeirMoM = FALSE; gridLength=11; gridRange=c(-5, 5); trend = c("none", "commonDispersion", "trendedDispersion")[1]; priorDf=10; span=0.3; mcCores=1; verbose=FALSE


dmEstimateTagwiseDisp <- function(dge, group = NULL, adjust = TRUE, mode = c("constrOptim", "constrOptim2", "constrOptim2G", "optim2", "optim2NM", "FisherScoring")[3], epsilon = 1e-05, maxIte = 1000, modeDisp = c("optimize", "optim", "constrOptim", "grid")[2], interval = c(0, 1e+5), tol = 1e-00,  initDisp = 10, initWeirMoM = TRUE, gridLength = 15, gridRange = c(-7, 7), trend = c("none", "commonDispersion", "trendedDispersion")[1], priorDf = 10, span = 0.3, mcCores = 20, verbose = FALSE, plot = FALSE){
  
  
  y <- dge$counts
  genes <- names(y)
  ngenes <- length(y)
	
	
  ### calculate mean expression of genes 
  meanExpr <- unlist(mclapply(seq(ngenes), function(g){ sum(y[[g]]) / ncol(y[[g]]) },  mc.cores=mcCores))  
  names(meanExpr) <- genes
  dge$meanExpr <- meanExpr
	
	
  
  if(is.null(group)) group <- dge$samples$group
  group <- as.factor(group)
  ngroups <- nlevels(group)
  lgroups <- levels(group)
  nlibs <- length(group)
  
  igroups <- list()
  for(gr in 1:ngroups){
    # gr=2
    igroups[[lgroups[gr]]] <- which(group == lgroups[gr])
    
  }
 
  
  ### Find optimized dispersion
  switch(modeDisp, 
         
         optimize={
           
           outList <- mclapply(seq(ngenes), function(g){
             # g = 1
             #     print(g)
             
             ### return NA if gene has 1 exon or observations in one sample in group (anyway this gene would not be fitted by dmFit)
             if(is.null(dmAdjustedProfileLikTG(gamma0 = interval[1] + (1-(sqrt(5) - 1)/2)*(interval[2]-interval[1]) , y = y[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjust = adjust, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose)))
               return(NA) 
             
             out <- optimize(f = dmAdjustedProfileLikTG, interval = interval,
                             y = y[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjust = adjust, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose,
                             maximum = TRUE, tol = tol) 
             
             return(out$maximum)  
             
           }, mc.cores=mcCores)
           
           
           names(outList) <- genes  
           dge$tagwiseDispersion <- unlist(outList)
           
         },
         
         
         optim={
           
      outList <- mclapply(seq(ngenes), function(g){
             # g = 12
						 if(verbose)
							 cat("gene", genes[g], "\n")
						 
             ### return NA if gene has 1 exon or observations in one sample in group (anyway this gene would not be fitted by dmFit)
             if(is.null(dmAdjustedProfileLikTG(gamma0 = initDisp, y = y[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjust = adjust, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose)))
               return(NA) 
             
             
             if(initWeirMoM)
               initDisp <- weirMoM(data = y[[g]], se=FALSE)
             #              print(initDisp)
             
             
             
             try( out <- optim(par = initDisp, fn = dmAdjustedProfileLikTG, gr = NULL, 
                          y = y[[g]], ngroups = ngroups, lgroups = lgroups, igroups = igroups, adjust = adjust, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose,
                          method = "L-BFGS-B", lower = 0 + epsilon, upper = 1e+15, control = list(fnscale = -1, factr = tol)) , silent = TRUE)
             
             #              print(out$par)
             
             
             return(c(out=out$par, init=initDisp))  
             
           }, mc.cores=mcCores)
           
           
           out <- do.call(rbind, outList)
           # rownames(out) <- genes  
           
           dge$tagwiseDispersion <- out[,"out"]
					 names(dge$tagwiseDispersion) <- genes
           dge$initDispersion <- out[,"init"]
           names(dge$initDispersion) <- genes
                  
         },



         constrOptim={
                    
           outList <- mclapply(seq(ngenes), function(g){
             # g = 1
             #     print(g)
             
             ### return NA if gene has 1 exon or observations in one sample in group (anyway this gene would not be fitted by dmFit)
             if(is.null(dmAdjustedProfileLikTG(gamma0 = initDisp, y = y[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjust = adjust, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose)))
               return(NA) 
             
             ui <- 1
             ci <- 0 + epsilon
             
             if(initWeirMoM)
               initDisp <- weirMoM(data = y[[g]], se=FALSE)
             
             
             out <- constrOptim(theta = initDisp, dmAdjustedProfileLikTG, grad = NULL, method = "Nelder-Mead",
                                ui=ui, ci=ci, control=list(fnscale = -1, reltol = tol), 
                                y = y[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjust = adjust, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose )
             
             
             return(c(out=out$par, init=initDisp)) 
             
           }, mc.cores=mcCores )
           
           out <- do.call(rbind, outList)
           # rownames(out) <- genes  
           
           dge$tagwiseDispersion <- out[,"out"]
					 names(dge$tagwiseDispersion) <- genes
           dge$initDispersion <- out[,"init"]
           names(dge$initDispersion) <- genes
           
         },
      
      
          grid={
        
            ### genrate spline dispersion
            splinePts <- seq(from = gridRange[1], to = gridRange[2], length = gridLength)
            splineDisp <- dge$commonDispersion * 2^splinePts
            
            ### calculate the likelihood for each gene at the spline dispersion points
            
            loglik0L <- mclapply(seq(length(y)), function(g){
              # g = 1
              
              ll <- numeric(gridLength)
              
              for(i in seq(gridLength)){
                # i = 10 
                out <- dmAdjustedProfileLikTG(gamma0 = splineDisp[i], y = y[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjust = adjust, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose)
                if(is.null(out))
                  return(NULL)
                
                ll[i] <- out
                
              }
              
              return(ll)
              
            }, mc.cores = mcCores)
            
            names(loglik0L) <- genes
            
            loglik0 <- do.call(rbind, loglik0L)
            
						if(is.null(loglik0)){
	            dge$tagwiseDispersion <- rep(dge$commonDispersion, ngenes)
	            names(dge$tagwiseDispersion) <- genes
							cat("** Tagwise dispersion: ", head(dge$tagwiseDispersion), "... \n")
							return(dge)
						}
						
						
            ngenes2 <- nrow(loglik0)
            genes2 <- rownames(loglik0)
            loglik <- loglik0 
            
            if(trend != "none"){
  
              switch(trend, 
                     commonDispersion={

                       moderation <- matrix(colMeans(loglik0), ngenes2, gridLength, byrow=TRUE)
                       
                     },
                     
                     trendedDispersion={
                       
                       o <- order(dge$meanExpr[genes2])
                       oo <- order(o)
                       width <- floor(span * ngenes2)
                       
                       moderation <- edgeR::movingAverageByCol(loglik0[o,], width=width)[oo,]
                       
                     })
              
              rownames(moderation) <- rownames(loglik0)
              
              priorN <- priorDf/(nlibs - ngroups) ### analogy to edgeR
              
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
            
            names(out) <- genes2
            
if(plot){
  dge$plotOutDisp <- dge$commonDispersion * 2^out
}
#### set common dispersion for genes that tagwise disp could not be calculated 
            dge$tagwiseDispersion <- rep(dge$commonDispersion, ngenes)
            names(dge$tagwiseDispersion) <- genes
            dge$tagwiseDispersion[genes2] <- dge$commonDispersion * 2^out
  
        
      })
  
  
  cat("** Tagwise dispersion: ", head(dge$tagwiseDispersion), "... \n")
  
  return(dge)
  
}






