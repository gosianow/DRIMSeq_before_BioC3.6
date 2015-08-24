#' @include class_show_utils.R
NULL

##############################################################

#' MatrixList container
#' 
#' @slot unlistData Matrix data
#' @slot partitioning List with indexes of each submatrix
#' @slot metadata Matrix of additional information
setClass("MatrixList", 
  representation(unlistData = "matrix", 
    partitioning = "list", 
    metadata = "matrix"))



##############################################################

MatrixList <- function(..., metadata){
  
  listData <- list(...)
  
  if (length(listData) == 1L && is.list(listData[[1L]]))
  listData <- listData[[1L]]
  
  if (length(listData) == 0L) {
    return(new("MatrixList"))
    
    } else {
      
      if (!all(sapply(listData, is, "matrix")))
      stop("all elements in '...' must be matrices!")
      
      unlistData <- do.call(rbind, listData)
      
      w <- sapply(listData, nrow)
      
      partitioning <- vector("list", length(w))
      
      partitioning[w != 0] <- split(1:nrow(unlistData), rep(1:length(w), w))
      
      if(!is.null(names(listData)))
      names(partitioning) <- names(listData)
      
      if(!missing(metadata))
      return(new("MatrixList", unlistData = unlistData, partitioning = partitioning, metadata = metadata))
      else
      return(new("MatrixList", unlistData = unlistData, partitioning = partitioning))
        
  }
    
}


##############################################################

#' @export
setMethod("length", "MatrixList", function(x){
  
  length(x@partitioning)
  
  })

#' @export
setMethod("names", "MatrixList", function(x){
  
  names(x@partitioning)
  
  })

##############################################################


setMethod("show", "MatrixList", function(object){
  
  nhead <- 2
  
  nl <- length(object)  
  cat(mode(object@unlistData),"MatrixList of length", nl,"\n")
  
  if(nl > 0){
    np <- min(nl, nhead)
    
    object <- object[1:np]
    
    if(is.null(names(object)))
      print_names <- paste0("[[", 1:np, "]]\n")
    else 
      print_names <- paste0("$", names(object), "\n")
    
    for(i in 1:np){
      
      cat(print_names[i])
      show_matrix(object[[i]])
      cat("\n")
      
    }
    
    if(np < nl){
      
      if(is.null(names(object)))
        cat(paste0("[[...]]\n"))
      else 
        cat(paste0("$...\n"))
      
    }
    
  }
  
  if(nrow(object@metadata) != 0){
    cat("\nwith metadata slot\n")
    show_matrix(object@metadata)
  }

  
})


##############################################################

#' @export
setMethod("[[", "MatrixList", function(x, i, ...){
  
  x@unlistData[x@partitioning[[i]], , drop = FALSE]
  
})


#' @export
setMethod("$", "MatrixList", function(x, name){
  
  x[[name]]
  
})


#' @export
setMethod("nrow", "MatrixList", function(x){
  
  nrow(x@unlistData)
  
})

#' @export
setMethod("ncol", "MatrixList", function(x){
  
  ncol(x@unlistData)
  
})

#' @export
setMethod("dim", "MatrixList", function(x){
  
  dim(x@unlistData)
  
})

##############################################################

#' @export
setMethod("width", "MatrixList", function(x){
  
  sapply(x@partitioning, length)
  
})


##############################################################

#' @export
setMethod("[", "MatrixList", function(x, i, j, ..., drop = TRUE){

  if(!missing(i)){
    
    if(!missing(j)){
      
      if(nrow(x@metadata) != 0)
      return(new("MatrixList", unlistData = x@unlistData[unlist(x@partitioning[i]), j, drop = FALSE], partitioning = relist(1:nrow(x@unlistData), x@partitioning[i]), metadata = x@metadata[i, , drop = FALSE]))
      else
      return(new("MatrixList", unlistData = x@unlistData[unlist(x@partitioning[i]), j, drop = FALSE], partitioning = relist(1:nrow(x@unlistData), x@partitioning[i]), metadata = x@metadata))
      

    }else{
      
      if(nrow(x@metadata) != 0)
      return(new("MatrixList", unlistData = x@unlistData[unlist(x@partitioning[i]), , drop = FALSE], partitioning = relist(1:nrow(x@unlistData), x@partitioning[i]), metadata = x@metadata[i, , drop = FALSE]))
      else
      return(new("MatrixList", unlistData = x@unlistData[unlist(x@partitioning[i]), , drop = FALSE], partitioning = relist(1:nrow(x@unlistData), x@partitioning[i]), metadata = x@metadata))   

    } 
    
  }else{
    
    if(!missing(j)){
      
      return(new("MatrixList", unlistData = x@unlistData[, j, drop = FALSE], partitioning = x@partitioning, metadata = x@metadata))

    }else{
      
      return(x)
      
    } 
    
  }

})







