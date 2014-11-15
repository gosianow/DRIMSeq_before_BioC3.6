##############################################################################
# adjustements to profile likelihood
##############################################################################


dmSQTLAdj <- function(gamma0, dgeSQTL, mcCores=20){
  
  y <- dgeSQTL$counts
  
  adj <- unlist(mclapply(seq(nrow(dgeSQTL$SNPs)), function(snp){  
    # snp = 1
    if(is.null(dgeSQTL$fit[[snp]])) return(NULL)
    
    g <- dgeSQTL$SNPs[snp, "gene_id"]
    gt <- dgeSQTL$genotypes[snp,]
    
    a <- dmSQTLAdjCROneGeneManyGroups(y = y[[g]][, !is.na(gt)], group = na.omit(gt), gamma0 = gamma0, piH = dgeSQTL$fit[[snp]]$piH) 
    
    return(a)
    
  }, mc.cores=mcCores))
  
  return(adj)
  
}


