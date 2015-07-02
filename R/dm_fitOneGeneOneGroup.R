##############################################################################
## estimate pi for given dispersion
##############################################################################


dm_fitOneGeneOneGroup <- function(y, gamma0, prop_mode = c("constrOptim", "constrOptimG", "FisherScoring")[2], prop_tol = 1e-12, verbose = FALSE){
  ### y must be exons vs. samples
	
  # NULL for filtered genes or genes with one exon
  if(any(dim(y) <= 1)) return(NULL)
  
  ### check for 0s in rows (exons)
  keep_row <- rowSums(y) > 0
  if(sum(keep_row) < 2) return(NULL) ## must be at least two exons or transcripts
  y <- y[keep_row, , drop=FALSE]
  
  ### check for 0s in cols (samples)
  keep_col <- colSums(y) > 0
  if(sum(keep_col) < 2) return(NULL) ## must be at least two samples in a condition
  y <- y[, keep_col, drop=FALSE]
  
  pi_init <- rowSums(y)/sum(y)
  k <- length(pi_init) ## k - number of exons
  N <- ncol(y) ## N - number of samples
	
  
  switch(prop_mode, 
         
         ### must have constraint for SUM pi = 1 --> sum(pi) < 1 + eps & sum(pi) > 1 - eps
         constrOptim = { ## for k-1 parameters
           # if(verbose) cat("\n gene:", colnames(y)[1], "gamma0:", gamma0, fill = T)
           
           ui <- rbind(diag(rep(1, k-1), k-1), diag(rep(-1, k-1), k-1), rep(-1, k-1))
           ci <- c(rep(0, k-1), rep(-1, k-1), -1 + .Machine$double.eps) 
           #         ui <- rbind(diag(rep(1, k-1)), diag(rep(-1, k-1)))
           #         ci <- c(rep(0, k-1), rep(-1, k-1))
           
           co <- constrOptim(pi_init[-k], f = dm_lik, grad = dm_score, ui = ui, ci = ci, control = list(fnscale = -1, reltol = prop_tol), gamma0 = gamma0, y = y)
           
           pi <- co$par
           lik <- co$value
           
           pi <- c(pi, 1-sum(pi))

         }, 
                 
         constrOptimG={ ## for k-1 parameters with Gamma functions
           # if(verbose) cat("\n gene:", colnames(y)[1], "gamma0:", gamma0, fill = T)
           
           ui <- rbind(diag(rep(1, k-1), k-1), diag(rep(-1, k-1), k-1), rep(-1, k-1))
           ci <- c(rep(0, k-1), rep(-1, k-1), -1 + .Machine$double.eps) 
           #         ui <- rbind(diag(rep(1, k-1)), diag(rep(-1, k-1)))
           #         ci <- c(rep(0, k-1), rep(-1, k-1))
           
           co <- constrOptim(pi_init[-k], f = dm_likG, grad = dm_scoreG, ui = ui, ci = ci, control = list(fnscale = -1, reltol = prop_tol), gamma0 = gamma0, y = y)
           
           pi <- co$par
           lik <- co$value
           
           pi <- c(pi, 1-sum(pi))

         }, 
         
         FisherScoring={
					 ### This mode is completely unprepared!!!
					 
					 plot = FALSE
					 epsilon = 1e-05
					 maxIte = 1000
           pi_initOrg <- pi_init <- rowSums(y)/sum(y)
           piMAX <- pi <- pi_init[-k]
           lik1 <- 0
           likMAX <- lik2 <- dmLogLikkm1(pi, gamma0, y)
           ite <- 1
           conv <- TRUE 
           
           if(plot){
             piX <- seq(0, 1, by = 0.01)
             piX <- piX[c(-1, -length(piX))]
             loglikY <- rep(0, length(piX))
             for(i in 1:length(loglikY))
               loglikY[i] <- dmLogLikkm1(piX[i], gamma0, y)
             
             plot(piX, loglikY, type="l", col="deeppink", lwd=4, main = gamma0)
             abline(v = pi, lty =  3)
             plot(piX, loglikY, type="l", col="deeppink", lwd=6, main = gamma0, xlim=c(pi - 0.3, pi + 0.3))
             abline(v = pi, lty =  3)
             points(pi, lik2, pch="0")
           }
           
           # Iterations
           while(conv & ite <= maxIte){
             if(verbose) cat("\n gene:", colnames(y)[1], "gamma0:", gamma0, fill = T)
             if(verbose) cat("ite:", ite, fill = T)
             if(verbose) cat("lik2-lik1:", lik2 - lik1, fill = T)
             
             if(abs(lik2 - lik1) < epsilon) conv <- FALSE
             
             scoreFun <- dm_scoreG(pi, gamma0, y)
             if(verbose) cat("scoreFun:", scoreFun, fill = T)
             
             # if(mode=="exp") invFIM <- dmInvExpFIMkm1(pi, gamma0, y)
             # else if(mode=="obs") invFIM <- dmInvObsFIMkm1(pi, gamma0, y)
						 invFIM <- dm_FIMobs(pi, gamma0, y)
						 
             if(verbose) cat("invFIM:", invFIM, fill = T)
             
             # Updates parameter estimates
             pi <- pi + invFIM %*% scoreFun
             
             ## check if pi is negative, then restart with new init params
             if(any(c(pi, 1-sum(pi)) <= 0 )){
               cat("* Negative pi:",c(pi, 1-sum(pi)), fill = T)
               ### generate new starting params
               #         randInit <- runif(k)
               #         pi_init <- randInit/sum(randInit)
               pi_init <- rdirichlet(1, pi_initOrg*10)
               while(any(pi_init < 1e-10))
                 pi_init <- rdirichlet(1, pi_initOrg*10)
               cat("pi_init:", pi_init, fill = T)
               pi <- pi_init[-k]
             }
             
             if(verbose) cat("pi:",c(pi, 1-sum(pi)), fill = T)
             
             lik1 <- lik2
             lik2 <- dm_likG(pi, gamma0, y)
             
             if(lik2 > likMAX){
               likMAX  <- lik2
               piMAX <- pi
             }
                         
             if(plot)
               points(pi, lik2, pch="|")
             if(verbose) cat("lik2:", lik2, fill = T)
             ite <- ite + 1
           }         
           
           if(ite > maxIte){
             pi <- piMAX
             lik1 <- likMAX
           }
           
           if(plot){
             points(pi, lik1, pch="*", cex=2, col="darkturquoise")
             dev.off()
           }
           
           pi <- c(pi, 1-sum(pi))
					 
				   #   ## check if the sum of pi equals 1
				   #   if(sum(pi) != 1){
				   #     cat("gene:", colnames(y)[1], "gamma0:", gamma0, fill = T)
				   #     cat("**** Not 1. sum(pi) = ", sum(pi), fill = T)
				   #     cat("**** pi:",pi, fill = T)
				   #   }
  
				   ## check if pi is negative
				   if(any(pi < 0 | pi >1 )){
				     cat("**** Negative pi:",pi, fill = T)
				   }
					 
         }) ## End of FisherScoring
  
  
  keep_row[keep_row] <- pi
  pi <- keep_row

  return(list(pi = pi, gamma0 = gamma0, lik = lik, df = k-1))
  
  
}







