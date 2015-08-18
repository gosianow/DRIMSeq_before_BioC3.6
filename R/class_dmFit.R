##############################################################
#' Object that contains Dirichlet-multinomial fitting results.
#' 
#' @slot proportions \code{\linkS4class{MatrixList}} of fitted proportions. Length of this list is equal to the number of genes. Each element is a matrix with rows corresponding to the features of a gene and columns to the conditions/grouping of samples.
#' @slot statistics DataFrame with per gene loglikelihoods \code{lik} and degrees of freedom \code{df} 
#' @importClassesFrom S4Vectors DataFrame
setClass("dmFit", 
  representation(
    proportions = "MatrixList",
    statistics = "DataFrame"))


setMethod("show", "dmFit", function(object){
  
  cat("An object of class", class(object), "\n")
  
  cat("\nSlot \"proportions\":\n")
  print(object@proportions)
  
  cat("\nSlot \"statistics\":\n")
  print(object@statistics)
  
  })
