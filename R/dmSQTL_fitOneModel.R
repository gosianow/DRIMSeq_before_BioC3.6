##############################################################################
# multiple group fitting 
##############################################################################

# counts=x@counts; genotypes=x@genotypes; dispersion = 38196.6; model = c("full", "null")[1]; prop_mode=c("constrOptim", "constrOptimG", "FisherScoring")[2]; prop_tol = 1e-12; verbose=FALSE; BPPARAM = BiocParallel::MulticoreParam(workers = 10)

# it returns a list of MatrixLists; unlistData = pi, metadata = stats


dmSQTL_fitOneModel <- function(counts, genotypes, dispersion, model = c("full", "null")[1], prop_mode=c("constrOptim", "constrOptimG", "FisherScoring")[2], prop_tol = 1e-12, verbose=FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  inds <- 1:length(counts)
  
  if(class(dispersion) == "numeric"){ 
    gamma0 <- relist(rep(dispersion, nrow(genotypes)), genotypes@partitioning)
  } else {
    gamma0 <- dispersion
  }
  
  lgroups_g <- c("0", "1", "2")
  ngroups_g <- 3
  
  switch(model, 
         
         full={
           
           if(verbose) cat("* Fitting full model.. \n")
           
           time <- system.time(fff <- BiocParallel::bplapply(inds, function(g){
             # g = 662
             
             y = counts[[g]]
             n_y <- nrow(y)
             snps = genotypes[[g]]
             n_snps <- nrow(snps)
             names_snps <- rownames(snps)
             
             pi <- matrix(NA, nrow = n_y * n_snps, ncol = ngroups_g, dimnames = list(rep(rownames(y), n_snps), lgroups_g))
             stats <- matrix(NA, n_snps, ngroups_g, dimnames = list(names_snps, lgroups_g))
             
             for(i in 1:n_snps){          
              # i = 29
              print(i)

              NAs <- is.na(snps[i, ]) | is.na(y[1, ])            
              yg <- y[, !NAs]             
              group <- snps[i, !NAs]
              group <- factor(group)
              ngroups <- nlevels(group)
              lgroups <- levels(group)
              nlibs <- length(group)
              
              igroups <- lapply(lgroups, function(gr){which(group == gr)})
              names(igroups) <- lgroups
              
              f <- dm_fitOneGeneManyGroups(y = yg, ngroups = ngroups, lgroups = lgroups, igroups = igroups, gamma0 = gamma0[[g]][i], prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)  
              
              ipi <- (i-1)*n_y + 1
              
              pi[ipi:(ipi+n_y-1), lgroups] <- f$pi
              stats[i, lgroups] <- f$stats
              
             }
             
             partitioning <- split(1:nrow(pi), factor(rep(1:n_snps, each = n_y)))
             names(partitioning) <- names_snps
             
             ff <- new("MatrixList", unlistData = pi, partitioning = partitioning, metadata = stats)
             
             
             return(ff)
              
           }, BPPARAM = BPPARAM))
           
           if(verbose) cat("Took ", time["elapsed"], " seconds.\n")
           names(fff) <- names(counts)  
           
           return(fff)
           
         }, 
         
         null={
           
           if(verbose) cat("* Fitting null model.. \n")
           
           time <- system.time(fff <- BiocParallel::bplapply(inds, function(g){
             # g = 1; y = counts[[g]]; snps = genotypes[[g]]
             
             y = counts[[g]]
             n_y <- nrow(y)
             snps = genotypes[[g]]
             n_snps <- nrow(snps)
             names_snps <- rownames(snps)
             
             pi <- matrix(NA, nrow = n_y * n_snps, ncol = 1, dimnames = list(rep(rownames(y), n_snps), "null"))
             stats <- matrix(NA, n_snps, 2, dimnames = list(names_snps, c("lik", "df")))
             
             
             for(i in 1:n_snps){          
              # i = 1
              # print(i)
              
              NAs <- is.na(snps[i, ]) | is.na(y[1, ])            
              yg <- y[, !NAs]
                           
              f <- dm_fitOneGeneOneGroup(y = yg, gamma0 = gamma0[[g]][i], prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)
              
              ipi <- (i-1)*n_y + 1
              
              pi[ipi:(ipi+n_y-1), ] <- f$pi
              stats[i, ] <- f$stats
              
             }
             
             partitioning <- split(1:nrow(pi), factor(rep(1:n_snps, each = n_y)))
             names(partitioning) <- names_snps
             
             ff <- new("MatrixList", unlistData = pi, partitioning = partitioning, metadata = stats)
             
             return(ff)
             
             
           }, BPPARAM = BPPARAM))
           
           if(verbose) cat("Took ", time["elapsed"], " seconds.\n")
           names(fff) <- names(counts)
           
           
           return(fff)
           
         })

}

