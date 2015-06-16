##############################################################################
# multiple group fitting 
##############################################################################


dmFit <- function(dge, model = c("full", "null")[1], dispersion = c("commonDispersion", "tagwiseDispersion")[1], modeProp=c("constrOptim2", "constrOptim2G", "FisherScoring")[2], tolProp = 1e-12, verbose=FALSE, BPPARAM = MulticoreParam(workers=1)){
  
	if(model == "full")
		cat("Fitting full model.. \n")
	
  geneList <-  names(dge$counts)
  
  
  ### dispersion / TODO: add errors in length of dispersion is not equal to length(genes)
  
  if(is.character(dispersion)){
    switch(dispersion, 
           commonDispersion = { gamma0 <- rep(dge$commonDispersion, length(geneList)); names(gamma0) <- geneList},
           tagwiseDispersion = { gamma0 <- dge$tagwiseDispersion } )
    
  } else {
    
    if(length(dispersion == 1)){
     gamma0 <- rep(dispersion, length(geneList)); names(gamma0) <- geneList
    } else {
     gamma0 <- dispersion
    }
      
  }

  time <- system.time(switch(model, 
         full = {
         	
				  	group <- dge$samples$group
				   group <- as.factor(group)
				   ngroups <- nlevels(group)
				   lgroups <- levels(group)
  
				   igroups <- list()
				   for(gr in 1:ngroups){
				     # gr=2
				     igroups[[lgroups[gr]]] <- which(group == lgroups[gr])
    
				   }
  
				   fit <- bplapply(geneList, function(g){  
				     # g = "ENSG00000135778"
				     # cat("Gene:", genes[g], fill = TRUE)
    
				     if(is.na(gamma0[g]))
				       return(NULL)
						 
						 yg <- dge$counts[[g]]
		
				     f <- dmOneGeneManyGroups(y = yg, ngroups = ngroups, lgroups = lgroups, igroups = igroups, gamma0 = gamma0[g], modeProp = modeProp, tolProp = tolProp, verbose = verbose)
		
				     return(f)
				   }, BPPARAM = BPPARAM)
  
				   names(fit) <- geneList
  
				   dge$fitFull <- fit
  
	
         },
					null = {
						
			
				   group <- factor(rep("Null", length(dge$samples$group)))
				   ngroups <- 1
				   lgroups <- "Null"
				   igroups <- list(Null = 1:length(dge$samples$group))
  
	
				   fit <- bplapply(geneList, function(g){  
				     # g = "ENSG00000135778"
				     # cat("Gene:", genes[g], fill = TRUE)
    
				     if(is.na(gamma0[g]))
				       return(NULL)
						 
						 yg <- dge$counts[[g]]
		
				     f <- dmOneGeneManyGroups(y = yg, ngroups = ngroups, lgroups = lgroups, igroups = igroups, gamma0 = gamma0[g], modeProp = modeProp, tolProp = tolProp, verbose = verbose)
		
				     return(f)
				   }, BPPARAM = BPPARAM)
  
				   names(fit) <- geneList
  
				   dge$fitNull <- fit
						
					}))
 
 if(model == "full")
 cat("Took ", time["elapsed"], " seconds.\n")
 	
  return(dge)
  
}


