##############################################################################
# calculate tagwise dispersions 
# dmEstimateTagwiseDisp, dmSQTLEstimateTagwiseDisp
##############################################################################

# group=NULL; adjust = TRUE; mode = "constrOptim2G"; epsilon = 1e-05; maxIte = 1000; interval = c(0, 1e+5); tol = 1e-00; mcCores=20; verbose=FALSE; modeDisp=c("optimize", "optim", "constrOptim")[2]; initDisp = 10; initWeirMoM = TRUE


dmEstimateTagwiseDisp <- function(dge, group=NULL, adjust = FALSE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=20, verbose=FALSE, modeDisp=c("optimize", "optim", "constrOptim")[1], initDisp = 10, initWeirMoM = FALSE){
  
  y <- dge$counts
  genes <- names(y)
  
  if(is.null(group)) group <- dge$samples$group
  group <- as.factor(group)
  ngroups <- nlevels(group)
  lgroups <- levels(group)
  
  igroups <- list()
  for(gr in 1:ngroups){
    # gr=2
    igroups[[lgroups[gr]]] <- which(group == lgroups[gr])
    
  }
  
  
  
  ### Find optimized dispersion
  switch(modeDisp, 
         
         optimize={
           
           
           outList <- mclapply(seq(length(y)), function(g){
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
           
           
           outList <- mclapply(seq(length(y)), function(g){
             # g = 12
#              print(g)
             

             ### return NA if gene has 1 exon or observations in one sample in group (anyway this gene would not be fitted by dmFit)
             if(is.null(dmAdjustedProfileLikTG(gamma0 = initDisp, y = y[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjust = adjust, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose)))
               return(NA) 
             
             
             if(initWeirMoM)
               initDisp <- weirMoM(data = y[[g]], se=FALSE)
             #              print(initDisp)
             
             
             
             try( out <- optim(par = initDisp, fn = dmAdjustedProfileLikTG, gr = NULL, 
                          y = y[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjust = adjust, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose,
                          method = "L-BFGS-B", lower = 0 + epsilon, upper = 1e+15, control=list(fnscale = -1)) , silent = TRUE)
             
             #              print(out$par)
             
             
             return(c(out=out$par, init=initDisp))  
             
           }, mc.cores=mcCores)
           
           
           out <- do.call(rbind, outList)
           rownames(out) <- genes  
           
           dge$tagwiseDispersion <- out[,"out"]
           dge$initDispersion <- out[,"init"]
           
           
           
         },
         constrOptim={
           
           
           outList <- mclapply(seq(length(y)), function(g){
             # g = 1
             #     print(g)
             
             ### return NA if gene has 1 exon or observations in one sample in group (anyway this gene would not be fitted by dmFit)
             if(is.null(dmAdjustedProfileLikTG(gamma0 = initDisp, y = y[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjust = adjust, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose)))
               return(NA) 
             
             ui <- 1
             ci <- 0 + epsilon
             
             if(initWeirMoM)
               initDisp <- weirMoM(data = y[[g]], se=FALSE)
             
             
             out <- constrOptim(theta = initDisp, dmAdjustedProfileLikTG, grad = NULL,
                                ui=ui, ci=ci, control=list(fnscale = -1), 
                                y = y[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjust = adjust, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose )
             
             
             return(c(out=out$par, init=initDisp)) 
             
           }, mc.cores=mcCores )
           
           out <- do.call(rbind, outList)
           rownames(out) <- genes  
           
           dge$tagwiseDispersion <- out[,"out"]
           dge$initDispersion <- out[,"init"]
           
           
         })
  
  
  
  cat("** Tagwise dispersion: ", head(dge$tagwiseDispersion), "... \n")
  
  return(dge)
  
}






