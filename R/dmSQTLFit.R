##############################################################################
# multiple group fitting 
##############################################################################


dmSQTLFit <- function(dgeSQTL, model = c("full", "null")[1], dispersion = c("commonDispersion", "tagwiseDispersion")[1], modeProp=c("constrOptim2", "constrOptim2G", "FisherScoring")[2], tolProp = 1e-12, verbose=FALSE, BPPARAM = MulticoreParam(workers=1)){
  
  if(is.character(dispersion)){
    switch(dispersion, 
           commonDispersion = { 
             gamma0 <- dgeSQTL$genotypes
             gamma0 <- lapply(gamma0, function(g){
               disp <- rep(dgeSQTL$commonDispersion, nrow(g))
               names(disp) <- rownames(g)
               return(disp)
             })
           },
           tagwiseDispersion = { gamma0 <- dgeSQTL$tagwiseDispersion } )
    
  } else {
    
    if(length(dispersion == 1)){
      gamma0 <- dgeSQTL$genotypes
      gamma0 <- lapply(gamma0, function(g){
        disp <- rep(dispersion, nrow(g))
        names(disp) <- rownames(g)
        return(disp)
      })
    } else {
      gamma0 <- dispersion
    }
    
  }
  
	geneList <- names(dgeSQTL$counts)
	
  switch(model, 
         full={
           
           fit <- bplapply(geneList, function(g){
             # g = geneList[1]; y = dgeSQTL$counts[[g]]; snps = dgeSQTL$genotypes[[g]]
             
             y = dgeSQTL$counts[[g]]
             snps = dgeSQTL$genotypes[[g]]
             # snps <- snps[1:min(nrow(snps), 5), , drop = FALSE]
             
             f <- list()
             
             for(i in 1:nrow(snps)){
               # i = 1
               if(is.na(gamma0[[g]][i]))
                 f[[i]] <- NULL
               
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
              
               
               f[[i]] <- dmOneGeneManyGroups(y = yg, ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
                                             gamma0 = gamma0[[g]][i], modeProp = modeProp, tolProp = tolProp, verbose = verbose)  
               
             }
             names(f) <- rownames(snps)
             
             return(f)
           }, BPPARAM = BPPARAM)
           
           names(fit) <- geneList
           dgeSQTL$fitFull <- fit
           
         }, 
         
         null={
           
           fit <- bplapply(geneList, function(g){
             # g = geneList[1]; y = dgeSQTL$counts[[g]]; snps = dgeSQTL$genotypes[[g]]
             
             y = dgeSQTL$counts[[g]]
             snps = dgeSQTL$genotypes[[g]]
             # snps <- snps[1:min(nrow(snps), 5), , drop = FALSE]
             
						 groupg <- factor(rep("Null", ncol(y)))
             ngroups <- 1
             lgroups <- "Null"	 
						 
             f <- list()
             
             for(i in 1:nrow(snps)){
               # i = 1
               if(is.na(gamma0[[g]][i]))
                 f[[i]] <- NULL
               
               NAs <- is.na(snps[i, ]) | is.na(y[1, ])            
               yg <- y[, !NAs]             
               group <- groupg[!NAs]
               nlibs <- sum(!NAs)              
               igroups <- list(Null = 1:nlibs)

               f[[i]] <- dmOneGeneManyGroups(y = yg, ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
                                             gamma0 = gamma0[[g]][i], modeProp = modeProp, tolProp = tolProp, verbose = verbose)
               
             }
             names(f) <- rownames(snps)
             
             return(f)
           }, BPPARAM = BPPARAM)
           
           names(fit) <- geneList
           
           dgeSQTL$fitNull <- fit
           
         })
  
  return(dgeSQTL)
  
}

