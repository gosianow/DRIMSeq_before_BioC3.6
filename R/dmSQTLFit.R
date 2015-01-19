##############################################################################
# multiple group fitting 
# dmFit, dmSQTLFit
##############################################################################




# model="full"; dispersion=3; mode="constrOptim2G"; epsilon = 1e-05; maxIte = 1000; verbose=FALSE; mcCores = 1


dmSQTLFit <- function(dgeSQTL, model="full", dispersion=c("commonDispersion", "tagwiseDispersion")[1], mode=c("constrOptim", "constrOptim2", "constrOptim2G", "optim2", "optim2NM", "FisherScoring")[3], epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 10){
  
  y <- dgeSQTL$counts
  
  if(is.character(dispersion)){
    switch(dispersion, 
           commonDispersion = { dgeSQTL$dispersion <- rep(dgeSQTL$commonDispersion, nrow(dgeSQTL$SNPs)) },
           tagwiseDispersion = { dgeSQTL$dispersion <- dgeSQTL$tagwiseDispersion } )
    
  } else {
    
    if(length(dispersion == 1)){
      dgeSQTL$dispersion <- rep(dispersion, nrow(dgeSQTL$SNPs))
    } else {
      dgeSQTL$dispersion <- dispersion
    }
    
  }
  
  gamma0 <- dgeSQTL$dispersion
  
  
  switch(model, 
         full={
                   
           fit <- mclapply(seq(nrow(dgeSQTL$SNPs)), function(snp){  
             # snp = 1
             # cat("SNP:", dgeSQTL$SNPs$SNP_id[snp], fill = TRUE)

						 if(is.na(gamma0[snp]))
							 return(NULL)

             NAs <- !is.na(dgeSQTL$genotypes[snp,]) & !is.na(y[[dgeSQTL$SNPs[snp, "gene_id"]]][1, ])
             
             y.g <- y[[dgeSQTL$SNPs[snp, "gene_id"]]][, NAs]
             
             group <- dgeSQTL$genotypes[snp, NAs]
             group <- as.factor(group)
             ngroups <- nlevels(group)
             lgroups <- levels(group)
             
             igroups <- list()
             for(gr in 1:ngroups){
               # gr=2
               igroups[[lgroups[gr]]] <- which(group == lgroups[gr])
               
             }
             
          
             
             f <- dmOneGeneManyGroups(y = y.g, ngroups = ngroups, lgroups = lgroups, igroups = igroups, gamma0 = gamma0[snp], mode = mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose)  
             
             
             return(f)
           }, mc.cores=mcCores)
           
           names(fit) <- dgeSQTL$SNPs$SNP_id
           
         }, 

         null={
           
           fit <- mclapply(seq(nrow(dgeSQTL$SNPs)), function(snp){  
						 
             # cat("SNP:", dgeSQTL$SNPs$SNP_id[snp], fill = TRUE)
						 
						 if(is.na(gamma0[snp]))
							 return(NULL)
						 
             NAs <- !is.na(dgeSQTL$genotypes[snp,]) & !is.na(y[[dgeSQTL$SNPs[snp, "gene_id"]]][1, ])
             
             y.g <- y[[dgeSQTL$SNPs[snp, "gene_id"]]][, NAs]
             
             f <- dmOneGeneGroup(y = y.g, gamma0 = gamma0[snp], mode = mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose)  

             return(f)
           }, mc.cores=mcCores)
           
           names(fit) <- dgeSQTL$SNPs$SNP_id
           
         })
  
  dgeSQTL$fit <- fit

  return(dgeSQTL)
  
}

