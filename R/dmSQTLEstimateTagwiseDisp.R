##############################################################################
# calculate tagwise dispersions 
##############################################################################

# group=NULL; adjust = FALSE; mode = "constrOptim2G"; epsilon = 1e-05; maxIte = 1000; modeDisp=c("optimize", "optim", "constrOptim", "grid")[4]; interval = c(0, 1e+5); tol = 1e-00;  initDisp = 10; initWeirMoM = FALSE; gridLength=11; gridRange=c(-5, 5); trend = c("none", "commonDispersion", "trendedDispersion")[2]; priorDf=10; span=0.3; mcCores=20; verbose=FALSE


dmSQTLEstimateTagwiseDisp <- function(dgeSQTL, adjust = TRUE, mode = c("constrOptim", "constrOptim2", "constrOptim2G", "optim2", "optim2NM", "FisherScoring")[3], epsilon = 1e-05, maxIte = 1000, modeDisp = c("optimize", "optim", "constrOptim", "grid")[2], interval = c(0, 1e+5), tol = 1e-00,  initDisp = 10, initWeirMoM = TRUE, gridLength = 15, gridRange = c(-7, 7), trend = c("none", "commonDispersion", "trendedDispersion")[1], priorDf = 10, span = 0.3, mcCores = 20, verbose = FALSE, plot = FALSE){
  
  
  y <- dgeSQTL$counts
  genes <- names(y)
  ngenes <- length(y)
  
 
  ### calculate mean expression of genes 
  meanExpr <- unlist(mclapply(seq(ngenes), function(g){ sum(y[[g]][,!is.na(y[[g]][1,])]) /  sum(!is.na(y[[g]][1,])) },  mc.cores=mcCores)) 
  names(meanExpr) <- genes
  dgeSQTL$meanExpr <- meanExpr
  
	
	
  ### Find optimized dispersion
  switch(modeDisp, 
         
         optimize={
           
           outList <- mclapply(seq(nrow(dgeSQTL$SNPs)), function(snp){
             # snp = 1
             # print(snp)
				     NAs <- !is.na(dgeSQTL$genotypes[snp,]) & !is.na(y[[dgeSQTL$SNPs[snp, "gene_id"]]][1, ])
             
				             y.g <- y[[dgeSQTL$SNPs[snp, "gene_id"]]][, NAs]
             
				             group <- dgeSQTL$genotypes[snp, NAs]
					   group <- as.factor(group)
					   ngroups <- nlevels(group)
					   lgroups <- levels(group)
					   nlibs <- length(group)
  
					   igroups <- list()
					   for(gr in 1:ngroups){
					     # gr=2
					     igroups[[lgroups[gr]]] <- which(group == lgroups[gr])
    
					   }
						 
             ### return NA if gene has 1 exon or observations in one sample in group (anyway this gene would not be fitted by dmFit)
             if(is.null(dmAdjustedProfileLikTG(gamma0 = interval[1] + (1-(sqrt(5) - 1)/2)*(interval[2]-interval[1]) , y = y.g, ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjust = adjust, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose)))
               return(NA) 
             
             out <- optimize(f = dmAdjustedProfileLikTG, interval = interval,
                             y = y.g, ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjust = adjust, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose,
                             maximum = TRUE, tol = tol) 
             
             return(out$maximum)  
             
           }, mc.cores=mcCores)
           
           
           names(outList) <- genes  
           dgeSQTL$tagwiseDispersion <- unlist(outList)
           
         },
         
         
         optim={
           
      outList <- mclapply(seq(nrow(dgeSQTL$SNPs)), function(snp){
             # snp = 1
             # print(snp)
				     NAs <- !is.na(dgeSQTL$genotypes[snp,]) & !is.na(y[[dgeSQTL$SNPs[snp, "gene_id"]]][1, ])
             
				             y.g <- y[[dgeSQTL$SNPs[snp, "gene_id"]]][, NAs]
             
				             group <- dgeSQTL$genotypes[snp, NAs]
					   group <- as.factor(group)
					   ngroups <- nlevels(group)
					   lgroups <- levels(group)
					   nlibs <- length(group)
  
					   igroups <- list()
					   for(gr in 1:ngroups){
					     # gr=2
					     igroups[[lgroups[gr]]] <- which(group == lgroups[gr])
    
					   }
      
						 
             ### return NA if gene has 1 exon or observations in one sample in group (anyway this gene would not be fitted by dmFit)
             if(is.null(dmAdjustedProfileLikTG(gamma0 = initDisp, y = y.g, ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjust = adjust, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose)))
               return(NA) 
             
             
             if(initWeirMoM)
               initDisp <- weirMoM(data = y.g, se=FALSE)
             #              print(initDisp)
             
             
             
             try( out <- optim(par = initDisp, fn = dmAdjustedProfileLikTG, gr = NULL, 
                          y = y.g, ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjust = adjust, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose,
                          method = "L-BFGS-B", lower = 0 + epsilon, upper = 1e+15, control=list(fnscale = -1)) , silent = TRUE)
             
             #              print(out$par)
             
             
             return(c(out=out$par, init=initDisp))  
             
           }, mc.cores=mcCores)
           
           
           out <- do.call(rbind, outList)
           rownames(out) <- genes  
           
           dgeSQTL$tagwiseDispersion <- out[,"out"]
           dgeSQTL$initDispersion <- out[,"init"]
                  
         },



         constrOptim={
                    
           outList <- mclapply(seq(nrow(dgeSQTL$SNPs)), function(snp){
             # snp = 1
             # print(snp)
				     NAs <- !is.na(dgeSQTL$genotypes[snp,]) & !is.na(y[[dgeSQTL$SNPs[snp, "gene_id"]]][1, ])
             
				             y.g <- y[[dgeSQTL$SNPs[snp, "gene_id"]]][, NAs]
             
				             group <- dgeSQTL$genotypes[snp, NAs]
					   group <- as.factor(group)
					   ngroups <- nlevels(group)
					   lgroups <- levels(group)
					   nlibs <- length(group)
  
					   igroups <- list()
					   for(gr in 1:ngroups){
					     # gr=2
					     igroups[[lgroups[gr]]] <- which(group == lgroups[gr])
    
					   }
             
             ### return NA if gene has 1 exon or observations in one sample in group (anyway this gene would not be fitted by dmFit)
             if(is.null(dmAdjustedProfileLikTG(gamma0 = initDisp, y = y.g, ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjust = adjust, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose)))
               return(NA) 
             
             ui <- 1
             ci <- 0 + epsilon
             
             if(initWeirMoM)
               initDisp <- weirMoM(data = y.g, se=FALSE)
             
             
             out <- constrOptim(theta = initDisp, dmAdjustedProfileLikTG, grad = NULL,
                                ui=ui, ci=ci, control=list(fnscale = -1), 
                                y = y.g, ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjust = adjust, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose )
             
             
             return(c(out=out$par, init=initDisp)) 
             
           }, mc.cores=mcCores )
           
           out <- do.call(rbind, outList)
           rownames(out) <- genes  
           
           dgeSQTL$tagwiseDispersion <- out[,"out"]
           dgeSQTL$initDispersion <- out[,"init"]
           
           
         },
      
      
          grid={
        
            ### genrate spline dispersion
            splinePts <- seq(from = gridRange[1], to = gridRange[2], length = gridLength)
            splineDisp <- dgeSQTL$commonDispersion * 2^splinePts
            
            ### calculate the likelihood for each gene at the spline dispersion points
            
            loglik0L <- mclapply(seq(nrow(dgeSQTL$SNPs)), function(snp){
             # snp = 1
             # print(snp)
						 
				     NAs <- !is.na(dgeSQTL$genotypes[snp,]) & !is.na(y[[dgeSQTL$SNPs[snp, "gene_id"]]][1, ])
             
				             y.g <- y[[dgeSQTL$SNPs[snp, "gene_id"]]][, NAs]
             
				             group <- dgeSQTL$genotypes[snp, NAs]
					   group <- as.factor(group)
					   ngroups <- nlevels(group)
					   lgroups <- levels(group)
					   nlibs <- length(group)
  
					   igroups <- list()
					   for(gr in 1:ngroups){
					     # gr=2
					     igroups[[lgroups[gr]]] <- which(group == lgroups[gr])
    
					   }
   		 
              ll <- numeric(gridLength)
              
              for(i in seq(gridLength)){
                # i = 1 
                out <- dmAdjustedProfileLikTG(gamma0 = splineDisp[i], y = y.g, ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjust = adjust, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose)
                if(is.null(out))
                  return(NULL)
                
                ll[i] <- out
                
              }
              
              return(ll)
              
            }, mc.cores = mcCores)
            
            names(loglik0L) <- dgeSQTL$SNPs$SNP_id
            
            loglik0 <- do.call(rbind, loglik0L)
            
            ngenes2 <- nrow(loglik0)
            genes2 <- rownames(loglik0)
            loglik <- loglik0 
            
            if(trend != "none"){
  
              switch(trend, 
                     commonDispersion={

                       moderation <- matrix(colMeans(loglik0), ngenes2, gridLength, byrow=TRUE)
                       
                     },
                     
                     trendedDispersion={
                       
                       o <- order(meanExpr[genes2])
                       oo <- order(o)
                       width <- floor(span * ngenes2)
                       
                       moderation <- edgeR::movingAverageByCol(loglik0[o,], width=width)[oo,]
                       
                     })
              
              rownames(moderation) <- rownames(loglik0)
              nlibs <- ncol(dgeSQTL$genotypes)
							ngroups <- 3
              priorN <- priorDf/(nlibs - ngroups) ### analogy to edgeR
              
              loglik <- loglik0 + priorN * moderation ### like in edgeR estimateTagwiseDisp
#               loglik <- (loglik0 + priorN * moderation)/(1 + priorN) ### like in edgeR dispCoxReidInterpolateTagwise

if(plot){
  
  dgeSQTL$plotSplineDisp <- splineDisp
  dgeSQTL$plotLoglik0 <- loglik0
  dgeSQTL$plotModeration <- moderation
  dgeSQTL$plotPriorN <- priorN
  dgeSQTL$plotLoglik <- loglik
  
}
              
            }
              

            out <- edgeR::maximizeInterpolant(splinePts, loglik)
              
            names(out) <- genes2
            
if(plot){
  dgeSQTL$plotOutDisp <- dgeSQTL$commonDispersion * 2^out
}

            dgeSQTL$tagwiseDispersion <- rep(NA, nrow(dgeSQTL$SNPs))
            names(dgeSQTL$tagwiseDispersion) <- dgeSQTL$SNPs$SNP_id
            dgeSQTL$tagwiseDispersion[genes2] <- dgeSQTL$commonDispersion * 2^out
  
        
      })
  
  
  cat("** Tagwise dispersion: ", head(dgeSQTL$tagwiseDispersion), "... \n")
  
  return(dgeSQTL)
  
}






