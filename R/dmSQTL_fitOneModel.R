##############################################################################
# multiple group fitting 
##############################################################################

# counts=x@counts; genotypes=x@genotypes; dispersion = 38196.6; model = c("full", "null")[1]; prop_mode=c("constrOptim", "constrOptimG", "FisherScoring")[2]; prop_tol = 1e-12; verbose=FALSE; BPPARAM = MulticoreParam(workers=10)

# it returns a list of list(pi - MatrixList, stats - matrix)


dmSQTL_fitOneModel <- function(counts, genotypes, dispersion, model = c("full", "null")[1], prop_mode=c("constrOptim", "constrOptimG", "FisherScoring")[2], prop_tol = 1e-12, verbose=FALSE, BPPARAM = MulticoreParam(workers=1)){
  
  gene_list <- names(counts)
  
  if(class(dispersion) == "numeric" && length(dispersion) == 1){ 
    gamma0 <- rep(dispersion, nrow(genotypes@unlistData))
    names(gamma0) <- rownames(genotypes@unlistData)
    # gamma0 <- IRanges::relist(gamma0, genotypes@partitioning) ### does not work, some problem with nchar...??
    gamma0 <- IRanges::relist(gamma0, IRanges::PartitioningByEnd(cumsum(IRanges::width(genotypes@partitioning))))
    names(gamma0) <- names(genotypes@partitioning)
  } else {
    gamma0 <- dispersion
  }
  
  
  switch(model, 
         
         full={
           
           cat("Fitting full model.. \n")
           
           time <- system.time(fff <- BiocParallel::bplapply(gene_list, function(g){
             # g = "ENSG00000131037.8"; y = counts[[g]]; snps = genotypes[[g]]
             
             y = counts[[g]]
             snps = genotypes[[g]]
             # snps <- snps[1:min(nrow(snps), 5), , drop = FALSE]
             
             ff <- lapply(rownames(snps), function(i){
               # i = rownames(snps)[3]
               
               NAs <- is.na(snps[i, ]) | is.na(y[1, ])            
               yg <- y[, !NAs]             
               group <- snps[i, !NAs]
               group <- factor(group)
               ngroups <- nlevels(group)
               lgroups <- levels(group)
               nlibs <- length(group)
               
               igroups <- lapply(lgroups, function(gr){which(group == gr)})
               names(igroups) <- lgroups
               
               f <- dm_fitOneGeneManyGroups(y = yg, ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
                                             gamma0 = gamma0[[g]][i], prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)  
               
               return(f)
               
             })
             
             names(ff) <- rownames(snps)
             
             pi <- MatrixList(lapply(ff, function(f){
              # print(f$pi)
              pp <- matrix(NA, nrow = nrow(f$pi), ncol = 3, dimnames = list(rownames(f$pi), 0:2)) 
              pp[, colnames(f$pi)] <- f$pi
              pp
              }))
             
             stats <- do.call(rbind, lapply(ff, function(f) f$stats))
             
             return(new("dmFit", proportions = pi, statistics = S4Vectors::DataFrame(stats, row.names = rownames(stats))))
             
             
           }, BPPARAM = BPPARAM))
           
           cat("Took ", time["elapsed"], " seconds.\n")
           names(fff) <- gene_list
           
           
           return(fff)
           
         }, 
         
         null={
           
           cat("Fitting null model.. \n")
           
           time <- system.time(fff <- BiocParallel::bplapply(gene_list, function(g){
             # g = gene_list[1]; y = counts[[g]]; snps = genotypes[[g]]
             
             y = counts[[g]]
             snps = genotypes[[g]]
             # snps <- snps[1:min(nrow(snps), 5), , drop = FALSE]
             
             groupg <- factor(rep("null", ncol(y)))
             ngroups <- 1
             lgroups <- "null"	 
             
             ff <- lapply(rownames(snps), function(i){
               # i = rownames(snps)[3]
               
               NAs <- is.na(snps[i, ]) | is.na(y[1, ])            
               yg <- y[, !NAs]             
               group <- groupg[!NAs]
               nlibs <- sum(!NAs)              
               igroups <- list(null = 1:nlibs)
               
               f <- dm_fitOneGeneManyGroups(y = yg, ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
                                         gamma0 = gamma0[[g]][i], prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)
               
               return(f)
               
             })
             
             names(ff) <- rownames(snps)
             
             pi <- MatrixList(lapply(ff, function(f) f$pi ))
             
             stats <- do.call(rbind, lapply(ff, function(f) f$stats ))
             
             return(new("dmFit", proportions = pi, statistics = S4Vectors::DataFrame(stats, row.names = rownames(stats))))
             
           }, BPPARAM = BPPARAM))
           
           cat("Took ", time["elapsed"], " seconds.\n")
           names(fff) <- gene_list
           
           
           return(fff)
           
         })

}

