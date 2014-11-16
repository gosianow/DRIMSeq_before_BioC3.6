
##############################################################################
# multiple group fitting 
# dmFit, dmSQTLFit
##############################################################################


# dge <- dgeDM; group=NULL; dispersion=NULL; mode="constrOptim2G"; epsilon = 1e-05; maxIte = 1000; verbose=FALSE; mcCores = 20


dmFit <- function(dge, group=NULL, dispersion=c("commonDispersion", "tagwiseDispersion")[2], mode=c("constrOptim", "constrOptim2", "constrOptim2G", "optim2", "optim2NM", "FisherScoring")[3], epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 20){
  
  y <- dge$counts
  genes <-  names(y)
  
  
  ### dispersion / TODO: add errors in length of dispersion is not equal to length(genes)
  
  if(is.character(dispersion)){
    switch(dispersion, 
           commonDispersion = { dge$dispersion <- rep(dge$commonDispersion, length(genes)) },
           tagwiseDispersion = { dge$dispersion <- dge$tagwiseDispersion } )
    
  } else {
    
    if(length(dispersion == 1)){
      dge$dispersion <- rep(dispersion, length(genes))
    } else {
      dge$dispersion <- dispersion
    }
      
  }

  gamma0 <- dge$dispersion
  
  if(is.null(group)) group <- dge$samples$group
  group <- as.factor(group)
  ngroups <- nlevels(group)
  lgroups <- levels(group)
  
  igroups <- list()
  for(gr in 1:ngroups){
    # gr=2
    igroups[[lgroups[gr]]] <- which(group == lgroups[gr])
    
  }
  
  fit <- mclapply(seq(length(y)), function(g){  
    # g = "ENSG00000135778"
    # cat("Gene:", genes[g], fill = TRUE)
    
    if(is.na(gamma0[g]))
      return(NULL)
    
    f <- dmOneGeneManyGroups(y[[g]], ngroups = ngroups, lgroups = lgroups, igroups = igroups, gamma0 = gamma0[g], mode = mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose)  
    
    return(f)
  }, mc.cores=mcCores)
  
  names(fit) <- genes
  
  
  dge$fit <- fit
  
  return(dge)
  
}


