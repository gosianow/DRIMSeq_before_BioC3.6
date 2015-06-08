##############################################################################
## estimate pi for given dispersion
##############################################################################


dmOneGeneGroup <- function(y, gamma0, modeProp = c("constrOptim2", "constrOptim2G", "FisherScoring")[2], tolProp = 1e-12, verbose = FALSE){
  ### y must be exons vs. samples
	
  # NULL for filtered genes or genes with one exon
  if(any(dim(y) <= 1)) return(NULL)
  
  ### check for 0s in rows (exons)
  keepRow <- rowSums(y) > 0
  if(sum(keepRow) < 2) return(NULL) ## must be at least two exons or transcripts
  y <- y[keepRow, , drop=FALSE]
  
  ### check for 0s in cols (samples)
  keepCol <- colSums(y) > 0
  if(sum(keepCol) < 2) return(NULL) ## must be at least two samples in a condition
  y <- y[, keepCol, drop=FALSE]
  
  piInit <- rowSums(y)/sum(y)
  k <- length(piInit) ## k - number of exons
  N <- ncol(y) ## N - number of samples
	
  
  switch(modeProp, 
         
         ### must have constraint for SUM pi = 1 --> sum(pi) < 1 + eps & sum(pi) > 1 - eps
         constrOptim2={ ## for k-1 parameters
           # if(verbose) cat("\n gene:", colnames(y)[1], "gamma0:", gamma0, fill = T)
           
           ui <- rbind(diag(rep(1, k-1), k-1), diag(rep(-1, k-1), k-1), rep(-1, k-1))
           ci <- c(rep(0, k-1), rep(-1, k-1), -1 + .Machine$double.eps) 
           #         ui <- rbind(diag(rep(1, k-1)), diag(rep(-1, k-1)))
           #         ci <- c(rep(0, k-1), rep(-1, k-1))
           
           co <- constrOptim(piInit[-k], f = dmLogLikkm1, grad = dmScoreFunkm1, ui = ui, ci = ci, control = list(fnscale = -1, reltol = tolProp), gamma0 = gamma0, y = y )
           
           piH <- co$par
           lik1 <- co$value
           
           # if(verbose) cat("piH:",c(piH, 1-sum(piH)), fill = T)
           # if(verbose) cat("lik1:", lik1, fill = T)
           
           piH <- c(piH, 1-sum(piH))
           # names(piH) <- names(piInit)
         }, 
         
         
         constrOptim2G={ ## for k-1 parameters with Gamma functions
           # if(verbose) cat("\n gene:", colnames(y)[1], "gamma0:", gamma0, fill = T)
           
           ui <- rbind(diag(rep(1, k-1), k-1), diag(rep(-1, k-1), k-1), rep(-1, k-1))
           ci <- c(rep(0, k-1), rep(-1, k-1), -1 + .Machine$double.eps) 
           #         ui <- rbind(diag(rep(1, k-1)), diag(rep(-1, k-1)))
           #         ci <- c(rep(0, k-1), rep(-1, k-1))
           
           co <- constrOptim(piInit[-k], f = dmLogLikGkm1, grad = dmScoreFunGkm1, ui = ui, ci = ci, control = list(fnscale = -1, reltol = tolProp), gamma0 = gamma0, y = y)
           
           piH <- co$par
           lik1 <- co$value
           
           # if(verbose) cat("piH:",c(piH, 1-sum(piH)), fill = T)
           # if(verbose) cat("lik1:", lik1, fill = T)
           
           piH <- c(piH, 1-sum(piH))
           # names(piH) <- names(piInit)
         }, 
         
         FisherScoring={
           # k-1 parameters for Fisher scoring
					 plot = FALSE
					 epsilon = 1e-05
					 maxIte = 1000
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
					 
         }) ## End of FisherScoring
  
  
  keepRow[keepRow] <- piH
  piH <- keepRow

  return(list(piH = piH, gamma0 = gamma0, logLik = lik1, df = k-1))
  
  
}







