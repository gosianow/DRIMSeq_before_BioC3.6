
##############################################################################
# estimates of pi -- Fisher scoring or constrOptim
# dmOneGeneGroup, dmOneGeneManyGroups, dmSQTLOneGeneManyGroups, 
##############################################################################
## estimate pi for given dispersion  with Fisher scoring -- one gene, one group

# y <- y[["g1"]]; gamma0 <- 1000;  mode = "obs"; epsilon=1e-05; maxIte = 1000; verbose = TRUE; plot = FALSE


dmOneGeneGroup <- function(y, gamma0, mode = c("constrOptim", "constrOptim2", "constrOptim2G", "optim2", "optim2NM", "FisherScoring")[3], epsilon = 1e-05, maxIte = 1000, verbose = FALSE, plot = FALSE){
  
  # NULL for filtered genes or genes with one exon
  if(dim(y)[1] <= 1) return(NULL)
  
  ### y must be exons vs. samples

  ### check for 0s in rows (exons)
  keepRow <- rowSums(y) > 0
  if(sum(keepRow) < 2) return(NULL) ## must be at least two exons or transcripts
  y <- y[keepRow, , drop=FALSE]
  df <- sum(keepRow)
  
  ### check for 0s in cols (samples)
  keepCol <- colSums(y) > 0
  if(sum(keepCol) < 2) return(NULL) ## must be at least two samples in a condition
  y <- y[, keepCol, drop=FALSE]
  
  piInit <- rowSums(y)/sum(y)
  k <- length(piInit) ## k - number of exons
  N <- ncol(y) ## N - number of samples
  
  if(plot){
    piX <- seq(0, 1, by = 0.01)
    piX <- piX[c(-1, -length(piX))]
    loglikY <- rep(0, length(piX))
    for(i in 1:length(loglikY))
      loglikY[i] <- dmLogLikkm1(piX[i], gamma0, y)
    
    plot(piX, loglikY, type="l", col="deeppink", lwd=4, main = gamma0)
    dev.off()
  }
  
  switch(mode, 
         
         constrOptim={ ## For k parameters

           ### must have constraint for SUM pi = 1 --> sum(pi) < 1 + eps & sum(pi) > 1 - eps
           ui <- rbind(diag(rep(1, k), k), diag(rep(-1, k), k), rep(1, k), rep(-1, k))
           ci <- c(rep(0, k), rep(-1, k), 1 - .Machine$double.eps, -1 - .Machine$double.eps)
           
           co <- constrOptim(piInit, f=dmLogLikG, grad=NULL, ui=ui, ci=ci, control = list(fnscale = -1), gamma0=gamma0, y=y )
           
           piH <- co$par
           lik1 <- co$value

         }, 
         
         ### must have constraint for SUM pi = 1 --> sum(pi) < 1 + eps & sum(pi) > 1 - eps
         constrOptim2={ ## for k-1 parameters
           # if(verbose) cat("\n gene:", colnames(y)[1], "gamma0:", gamma0, fill = T)
           
           ui <- rbind(diag(rep(1, k-1), k-1), diag(rep(-1, k-1), k-1), rep(-1, k-1))
           ci <- c(rep(0, k-1), rep(-1, k-1), -1 + .Machine$double.eps) # used to be with epsilon = 1e-5  but "initial value is not in the interior of the feasible region" because it is too big
           #         ui <- rbind(diag(rep(1, k-1)), diag(rep(-1, k-1)))
           #         ci <- c(rep(0, k-1), rep(-1, k-1))
           
           co <- constrOptim(piInit[-k], f=dmLogLikkm1, grad=dmScoreFunkm1, ui=ui, ci=ci, control=list(fnscale = -1), gamma0=gamma0, y=y )
#                co <- constrOptim(piInit[-k], f=dmLogLikkm1, grad=NULL, ui=ui, ci=ci, control=list(fnscale = -1), gamma0=gamma0, y=y )
           
           piH <- co$par
           lik1 <- co$value
           
           # if(verbose) cat("piH:",c(piH, 1-sum(piH)), fill = T)
           # if(verbose) cat("lik1:", lik1, fill = T)
           
           piH <- c(piH, 1-sum(piH))
           names(piH) <- names(piInit)
         }, 
         
         
         constrOptim2G={ ## for k-1 parameters with Gamma functions
           # if(verbose) cat("\n gene:", colnames(y)[1], "gamma0:", gamma0, fill = T)
           
           ui <- rbind(diag(rep(1, k-1), k-1), diag(rep(-1, k-1), k-1), rep(-1, k-1))
           ci <- c(rep(0, k-1), rep(-1, k-1), -1 + .Machine$double.eps) # used to be with epsilon = 1e-5  but "initial value is not in the interior of the feasible region" because it is too big
           #         ui <- rbind(diag(rep(1, k-1)), diag(rep(-1, k-1)))
           #         ci <- c(rep(0, k-1), rep(-1, k-1))
           
           co <- constrOptim(piInit[-k], f=dmLogLikGkm1, grad=dmScoreFunGkm1, ui=ui, ci=ci, control=list(fnscale = -1), gamma0=gamma0, y=y )
           
#            co <- constrOptim(piInit[-k], f=dmLogLikGkm1, grad=NULL, ui=ui, ci=ci, control=list(fnscale = -1), gamma0=gamma0, y=y )
           
           piH <- co$par
           lik1 <- co$value
           
           # if(verbose) cat("piH:",c(piH, 1-sum(piH)), fill = T)
           # if(verbose) cat("lik1:", lik1, fill = T)
           
           piH <- c(piH, 1-sum(piH))
           names(piH) <- names(piInit)
         }, 
         
         
         optim2={
           # if(verbose) cat("\n gene:", colnames(y)[1], "gamma0:", gamma0, fill = T)
           
           o <- optim(par = piInit[-k], fn = dmLogLikkm1, gr = dmScoreFunkm1, 
                      gamma0=gamma0, y=y,
                      method = "L-BFGS-B", lower = 0 + epsilon, upper = 1 - epsilon, control=list(fnscale = -1))
           
           piH <- o$par
           lik1 <- o$value
           piH <- c(piH, 1-sum(piH))
           names(piH) <- names(piInit)
           # if(verbose) cat("piH:", piH, fill = T)
         }, 
         
         
         optim2NM={
           # if(verbose) cat("\n gene:", colnames(y)[1], "gamma0:", gamma0, fill = T)
           o <- optim(par = piInit[-k], fn = dmLogLikkm1, gr = NULL, 
                      gamma0=gamma0, y=y,
                      method = c("Nelder-Mead", "BFGS", "CG", "SANN")[1], control=list(fnscale = -1))
           
           piH <- o$par
           lik1 <- o$value
           piH <- c(piH, 1-sum(piH))
           names(piH) <- names(piInit)
           # if(verbose) cat("piH:", piH, fill = T)
         }, 
         
         
         FisherScoring={
           # k-1 parameters for Fisher scoring
           piInitOrg <- piInit <- rowSums(y)/sum(y)
           piMAX <- piH <- piInit[-k]
           lik1 <- 0
           likMAX <- lik2 <- dmLogLikkm1(piH, gamma0, y)
           ite <- 1
           conv <- TRUE 
           
           if(plot){
             piX <- seq(0, 1, by = 0.01)
             piX <- piX[c(-1, -length(piX))]
             loglikY <- rep(0, length(piX))
             for(i in 1:length(loglikY))
               loglikY[i] <- dmLogLikkm1(piX[i], gamma0, y)
             
             plot(piX, loglikY, type="l", col="deeppink", lwd=4, main = gamma0)
             abline(v = piH, lty =  3)
             plot(piX, loglikY, type="l", col="deeppink", lwd=6, main = gamma0, xlim=c(piH - 0.3, piH + 0.3))
             abline(v = piH, lty =  3)
             points(piH, lik2, pch="0")
           }
           
           # Iterations
           while(conv & ite <= maxIte){
             if(verbose) cat("\n gene:", colnames(y)[1], "gamma0:", gamma0, fill = T)
             if(verbose) cat("ite:", ite, fill = T)
             if(verbose) cat("lik2-lik1:", lik2 - lik1, fill = T)
             
             if(abs(lik2 - lik1) < epsilon) conv <- FALSE
             
             scoreFun <- dmScoreFunkm1(piH, gamma0, y)
             if(verbose) cat("scoreFun:", scoreFun, fill = T)
             
             # if(mode=="exp") invFIM <- dmInvExpFIMkm1(piH, gamma0, y)
             # else if(mode=="obs") invFIM <- dmInvObsFIMkm1(piH, gamma0, y)
						 invFIM <- dmInvObsFIMkm1(piH, gamma0, y)
						 
             if(verbose) cat("invFIM:", invFIM, fill = T)
             
             # Updates parameter estimates
             piH <- piH + invFIM %*% scoreFun
             
             ## check if pi is negative, then restart with new init params
             if(any(c(piH, 1-sum(piH)) <= 0 )){
               cat("* Negative piH:",c(piH, 1-sum(piH)), fill = T)
               ### generate new starting params
               #         randInit <- runif(k)
               #         piInit <- randInit/sum(randInit)
               piInit <- rdirichlet(1, piInitOrg*10)
               while(any(piInit < 1e-10))
                 piInit <- rdirichlet(1, piInitOrg*10)
               cat("piInit:", piInit, fill = T)
               piH <- piInit[-k]
             }
             
             if(verbose) cat("piH:",c(piH, 1-sum(piH)), fill = T)
             
             lik1 <- lik2
             lik2 <- dmLogLikkm1(piH, gamma0, y)
             
             if(lik2 > likMAX){
               likMAX  <- lik2
               piMAX <- piH
             }
                         
             if(plot)
               points(piH, lik2, pch="|")
             if(verbose) cat("lik2:", lik2, fill = T)
             ite <- ite + 1
           }         
           
           if(ite > maxIte){
             piH <- piMAX
             lik1 <- likMAX
           }
           
           if(plot){
             points(piH, lik1, pch="*", cex=2, col="darkturquoise")
             dev.off()
           }
           
           piH <- c(piH, 1-sum(piH))
           names(piH) <- names(piInit) 
         }) ## End of FisherScoring
  
  
  keepRow[keepRow] <- piH
  piH <- keepRow
  
  #   ## check if the sum of pi equals 1
  #   if(sum(piH) != 1){
  #     cat("gene:", colnames(y)[1], "gamma0:", gamma0, fill = T)
  #     cat("**** Not 1. sum(piH) = ", sum(piH), fill = T)
  #     cat("**** piH:",piH, fill = T)
  #   }
  
  ## check if pi is negative
  if(any(piH < 0 | piH >1 )){
    cat("**** Negative piH:",piH, fill = T)
  }
  
  return(list(piH = as.matrix(piH), gamma0 = gamma0, logLik = lik1, df=df-1))
  
  
}







