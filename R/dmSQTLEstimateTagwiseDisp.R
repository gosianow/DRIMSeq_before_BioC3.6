##############################################################################
# calculate tagwise dispersions 
##############################################################################

# adjustDisp = TRUE; modeDisp = c("optimize", "optim", "constrOptim", "grid")[2]; intervalDisp = c(0, 1e+5); tolDisp = 1e-08;  initDisp = 10; initWeirMoMDisp = TRUE; gridLengthDisp = 15; gridRangeDisp = c(-7, 7); trendDisp = c("none", "commonDispersion", "trendedDispersion")[1]; priorDfDisp = 10; spanDisp = 0.3; modeProp = c( "constrOptim2", "constrOptim2G", "FisherScoring")[2]; tolProp = 1e-12; verbose = FALSE; plot = FALSE; BPPARAM = MulticoreParam(workers=1)

dmSQTLEstimateTagwiseDisp <- function(dgeSQTL, adjustDisp = TRUE, modeDisp = c("optimize", "optim", "constrOptim", "grid")[2], intervalDisp = c(0, 1e+5), tolDisp = 1e-08,  initDisp = 10, initWeirMoMDisp = TRUE, gridLengthDisp = 15, gridRangeDisp = c(-7, 7), trendDisp = c("none", "commonDispersion", "trendedDispersion")[1], priorDfDisp = 10, spanDisp = 0.3, modeProp = c( "constrOptim2", "constrOptim2G", "FisherScoring")[2], tolProp = 1e-12, verbose = FALSE, plot = FALSE, BPPARAM = MulticoreParam(workers=1)){
  
  
  ### calculate mean expression of genes 
  meanExpr <- unlist(bplapply(dgeSQTL$counts, function(g){ mean(colSums(g), na.rm = TRUE) }, BPPARAM = BPPARAM))	
  dgeSQTL$meanExpr <- meanExpr
  geneList <- names(dgeSQTL$counts)
  
  switch(modeDisp, 
         
         optimize={
           
           dispList <- bplapply(geneList, function(g){
             # g = geneList[1]; y = dgeSQTL$counts[[g]]; snps = dgeSQTL$genotypes[[g]]
             
             y = dgeSQTL$counts[[g]]
             snps = dgeSQTL$genotypes[[g]]
             # snps <- snps[1:min(nrow(snps), 5), , drop = FALSE]
             
             disp <- rep(0, nrow(snps))
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
               
               igroups <- list()
               for(gr in 1:ngroups){
                 # gr=2
                 igroups[[lgroups[gr]]] <- which(group == lgroups[gr])
                 
               }
               
               ### return NA if gene has 1 exon or observations in one sample in group (anyway this gene would not be fitted by dmFit)
               gamma0 <- intervalDisp[1] + (1-(sqrt(5) - 1)/2)*(intervalDisp[2]-intervalDisp[1])
               if(is.null(dmAdjustedProfileLikTG(gamma0 = gamma0, y = yg, ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjustDisp = adjustDisp, modeProp = modeProp, tolProp = tolProp, verbose = verbose))){
                 disp[i] <- NA
                 next
                 
               }
               
               
               optimum <- optimize(f = dmAdjustedProfileLikTG, interval = intervalDisp,
                                   y = yg, ngroups=ngroups, lgroups=lgroups, igroups=igroups, 
                                   adjustDisp = adjustDisp, modeProp = modeProp, tolProp = tolProp, verbose = verbose,
                                   maximum = TRUE, tol = tolDisp) 
               disp[i] <- optimum$maximum			 
               
             }
             
             return(disp)  
             
           }, BPPARAM = BPPARAM)
           
           names(dispList) <- geneList
           dgeSQTL$tagwiseDispersion <- dispList
           
         },
         
         
         
         
         optim={
           
           dispList <- bplapply(geneList, function(g){
             # g = geneList[1]; y = dgeSQTL$counts[[g]]; snps = dgeSQTL$genotypes[[g]]
             
             y = dgeSQTL$counts[[g]]
             snps = dgeSQTL$genotypes[[g]]
             # snps <- snps[1:min(nrow(snps), 5), , drop = FALSE]
             
             disp <- rep(0, nrow(snps))
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
               
               igroups <- list()
               for(gr in 1:ngroups){
                 # gr=2
                 igroups[[lgroups[gr]]] <- which(group == lgroups[gr])
                 
               }
               
               ### return NA if gene has 1 exon or observations in one sample in group (anyway this gene would not be fitted by dmFit)
               gamma0 <- initDisp
               if(is.null(dmAdjustedProfileLikTG(gamma0 = gamma0, y = yg, ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjustDisp = adjustDisp, modeProp = modeProp, tolProp = tolProp, verbose = verbose))){
                 disp[i] <- NA
                 next
               }
               
               if(initWeirMoMDisp)
                 initDisp <- weirMoM(data = yg, se=FALSE)
               
               try( optimum <- optim(par = initDisp, fn = dmAdjustedProfileLikTG, gr = NULL, 
                                     y = yg, ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
                                     adjustDisp = adjustDisp, modeProp = modeProp, tolProp = tolProp, verbose = verbose,
                                     method = "L-BFGS-B", lower = 1e-2, upper = 1e+10, control = list(fnscale = -1, factr = tolDisp)), silent = TRUE )
               
               disp[i] <- optimum$par			 
               
             }
             
             return(disp)  
             
           }, BPPARAM = BPPARAM)
           
           names(dispList) <- geneList
           dgeSQTL$tagwiseDispersion <- dispList
           
         },
         
         
         
         constrOptim={
           
           dispList <- bplapply(geneList, function(g){
             # g = geneList[1]; y = dgeSQTL$counts[[g]]; snps = dgeSQTL$genotypes[[g]]
             
             y = dgeSQTL$counts[[g]]
             snps = dgeSQTL$genotypes[[g]]
             # snps <- snps[1:min(nrow(snps), 5), , drop = FALSE]
             
             disp <- rep(0, nrow(snps))
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
               
               igroups <- list()
               for(gr in 1:ngroups){
                 # gr=2
                 igroups[[lgroups[gr]]] <- which(group == lgroups[gr])
                 
               }
               
               ### return NA if gene has 1 exon or observations in one sample in group (anyway this gene would not be fitted by dmFit)
               gamma0 <- initDisp
               if(is.null(dmAdjustedProfileLikTG(gamma0 = gamma0, y = yg, ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjustDisp = adjustDisp, modeProp = modeProp, tolProp = tolProp, verbose = verbose))){
                 disp[i] <- NA
                 next
                 
               }
               
               ui <- 1
               ci <- 1e-8
               
               if(initWeirMoMDisp)
                 initDisp <- weirMoM(data = yg, se=FALSE)
               
               
               optimum <- constrOptim(theta = initDisp, dmAdjustedProfileLikTG, grad = NULL, method = "Nelder-Mead",
                                      ui=ui, ci=ci, control = list(fnscale = -1, reltol = tolDisp), 
                                      y = yg, ngroups=ngroups, lgroups=lgroups, igroups=igroups, 
                                      adjustDisp = adjustDisp, modeProp = modeProp, tolProp = tolProp, verbose = verbose )
               
               disp[i] <- optimum$par			 
               
             }
             
             return(disp)  
             
           }, BPPARAM = BPPARAM)
           
           names(dispList) <- geneList
           dgeSQTL$tagwiseDispersion <- dispList
           
           
         },
         
         
         grid={
           
           ### genrate spline dispersion
           splinePts <- seq(from = gridRangeDisp[1], to = gridRangeDisp[2], length = gridLengthDisp)
           splineDisp <- dgeSQTL$commonDispersion * 2^splinePts
           
           ### calculate the likelihood for each gene at the spline dispersion points
           
           loglik0List <- bplapply(geneList, function(g){
             # g = geneList[1]; y = dgeSQTL$counts[[g]]; snps = dgeSQTL$genotypes[[g]]
             
             y = dgeSQTL$counts[[g]]
             snps = dgeSQTL$genotypes[[g]]
             # snps <- snps[1:min(nrow(snps), 5), , drop = FALSE]
             
             ll <- matrix(0, nrow(snps), gridLengthDisp)
             rownames(ll) <- paste0(g, "-",rownames(snps))
             
             for(i in 1:nrow(snps)){
               # i = 1
               
               NAs <- is.na(snps[i, ]) | is.na(y[1, ])            
               yg <- y[, !NAs]             
               group <- snps[i, !NAs]
               group <- factor(group)
               ngroups <- nlevels(group)
               lgroups <- levels(group)
               nlibs <- length(group)
               
               igroups <- list()
               for(gr in 1:ngroups){
                 # gr=2
                 igroups[[lgroups[gr]]] <- which(group == lgroups[gr])
                 
               }
               
               for(j in seq(gridLengthDisp)){
                 # j = 1 
                 out <- dmAdjustedProfileLikTG(gamma0 = splineDisp[j], y = yg, ngroups=ngroups, lgroups=lgroups, igroups=igroups, 
                                               adjustDisp = adjustDisp, modeProp = modeProp, tolProp = tolProp, verbose = verbose)
                 if(is.null(out)){
                   ll[i, ] <- NA
                   break
                 }
                 
                 ll[i, j] <- out
                 
               } # j
               
             } # i
             
             return(ll)
             
           }, BPPARAM = BPPARAM)
           
 
           loglik0 <- do.call(rbind, loglik0List)
           genesAll <- rownames(loglik0)
           loglik0 <- loglik0[complete.cases(loglik0), , drop = FALSE]
           genesComplete <- rownames(loglik0)
           loglik <- loglik0 
           
           if(trendDisp != "none"){
             
             switch(trendDisp, 
                    commonDispersion={
                      
                      moderation <- matrix(colMeans(loglik0), nrow(loglik0), gridLengthDisp, byrow=TRUE)
                      
                    },
                    
                    trendedDispersion={
                      
                      o <- order(meanExpr[limma::strsplit2(genesComplete, "-")[, 1]])
                      oo <- order(o)
                      width <- floor(spanDisp * nrow(loglik0))
                      
                      moderation <- edgeR::movingAverageByCol(loglik0[o,], width=width)[oo,]
                      
                    })
             
             rownames(moderation) <- genesComplete
             nlibs <- ncol(snps)
             ngroups <- 2
             priorN <- priorDfDisp/(nlibs - ngroups) ### analogy to edgeR
             
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
           names(out) <- genesComplete
           
           if(plot){
             dgeSQTL$plotOutDisp <- dgeSQTL$commonDispersion * 2^out
           }
           
           
           dispAll <- rep(NA, length(genesAll))
           names(dispAll) <- genesAll
           dispAll[genesComplete] <- dgeSQTL$commonDispersion * 2^out
           
           
           genesAll <- limma::strsplit2(genesAll, "-")
           names(dispAll) <- genesAll[, 2]
           
           dispAll <- split(dispAll, genesAll[, 1])
           dispAll <- dispAll[geneList]
           
           dgeSQTL$tagwiseDispersion <- dispAll
       
         })
  
  
  cat("** Tagwise dispersion done! \n")
  
  return(dgeSQTL)
  
}






