##############################################################################
# multiple group fitting 
# dmFit, dmSQTLFit
##############################################################################




# model="full"; dispersion=3000; mode="constrOptim2"; epsilon = 1e-05; maxIte = 1000; verbose=FALSE; mcCores = 20


dmSQTLFit <- function(dgeSQTL, model="full", dispersion=NULL, mode="constrOptim2", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 20){
  
  y <- dgeSQTL$counts
  
  if(is.null(dispersion)) dispersion <- dgeSQTL$commonDispersion
  gamma0 <- dgeSQTL$dispersion <- dispersion

  switch(model, 
         full={
                   
           fit <- mclapply(seq(nrow(dgeSQTL$SNPs)), function(snp){  
             # snp = 1
             # cat("SNP:", dgeSQTL$SNPs$SNP_id[snp], fill = TRUE)

             f <- dmSQTLOneGeneManyGroups(y = y[[dgeSQTL$SNPs[snp, "gene_id"]]][, !is.na(dgeSQTL$genotypes[snp,])], group = na.omit(dgeSQTL$genotypes[snp,]), gamma0 = gamma0, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose)  
#              gc()
             return(f)
           }, mc.cores=mcCores)
           
           names(fit) <- dgeSQTL$SNPs$SNP_id
           
         }, 
         null={
           
           fit <- mclapply(seq(nrow(dgeSQTL$SNPs)), function(snp){  
             f <- dmOneGeneGroup(y[[dgeSQTL$SNPs[snp, "gene_id"]]][, !is.na(dgeSQTL$genotypes[snp,])], gamma0 = gamma0, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose)  
#              gc()
             return(f)
           }, mc.cores=mcCores)
           
           names(fit) <- dgeSQTL$SNPs$SNP_id
           
         })
  
  dgeSQTL$fit <- fit

  return(dgeSQTL)
  
}

