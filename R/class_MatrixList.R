#' @include class_show_utils.R
NULL

##############################################################

#' MatrixList object
#' 
#' A MatrixList object is a 'list' container for matrices which have the same number of columns but varying number of rows. Columns of these matrices may or may not have the same interpretation. Additionally, one can store an extra information corresponding to each of the matrices in 'metadata' matrix. Below, you can find a list of methods that are available for MatrixList object. 
#' 
#' @slot unlistData A matrix which is a row binding of all the stored matrices.
#' @slot partitioning A list of indexes which defines the row partitioning of 'unlistData' matrix into the original submatrices.
#' @slot metadata A matrix of additional information where each row corresponds to one of the submatrices.
setClass("MatrixList", 
         representation(unlistData = "matrix", 
                        partitioning = "list", 
                        metadata = "matrix"))



setValidity("MatrixList", function(object){
  # has to return TRUE when valid object!
  
  partitioning_unlist <- unlist(object@partitioning)
  
  if(length(partitioning_unlist) == nrow(object@unlistData))
    out <- TRUE
  else
    return(paste0("Unequal lengths of partitioning indexes and rows in unlistData: ", length(partitioning_unlist), " and ", nrow(object@unlistData)))

  if(nrow(object@metadata) > 0){
    
    if(nrow(object@metadata) == length(object@partitioning))
      out <- TRUE
    else
      return(paste0("Unequal lengths of partitioning and metadata: ", length(object@partitioning), " and ", nrow(object@metadata)))
    
  }
  
  return(out)
  
})


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
    
    inds <- 1:nrow(unlistData)
    names(inds) <- rownames(unlistData)
    
    partitioning[w != 0] <- split(inds, rep(1:length(w), w))
    
    if(!is.null(names(listData)))
      names(partitioning) <- names(listData)
    
    if(!missing(metadata))
      return(new("MatrixList", unlistData = unlistData, partitioning = partitioning, metadata = metadata))
    else
      return(new("MatrixList", unlistData = unlistData, partitioning = partitioning))
    
  }
  
}


##############################################################


setMethod("show", "MatrixList", function(object){
  
  nhead <- 2
  
  nl <- length(object)  
  cat(mode(object@unlistData),"MatrixList of length", nl,"\n")
  
  if(nl > 0){
    np <- min(nl, nhead)
    
    object_sub <- object[1:np]
    
    if(is.null(names(object_sub)))
      print_names <- paste0("[[", 1:np, "]]\n")
    else 
      print_names <- paste0("$", names(object_sub), "\n")
    
    for(i in 1:np){
      # i = 1
      cat(print_names[i])
      show_matrix(object_sub[[i]])
      cat("\n")
      
    }
    
    if(np < nl){
      
      if(is.null(names(object_sub)))
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

#' @rdname MatrixList-class
#' @param x MatrixList object.
#' @export
setMethod("names", "MatrixList", function(x){
  
  names(x@partitioning)
  
})


#' @rdname MatrixList-class
#' @export
setMethod("names<-", "MatrixList", function(x, value){
  
  names(x@partitioning) <- value
  x
  
  # partitioning <- x@partitioning
  # names(partitioning) <- value
  # return(new("MatrixList", unlistData = x@unlistData, partitioning = partitioning, metadata = x@metadata))
  
})


#' @rdname MatrixList-class
#' @export
setMethod("rownames", "MatrixList", function(x){
  
  rownames(x@unlistData)
  
})


#' @rdname MatrixList-class
#' @param value Values used for the replacement of corresponding attributes.   
#' @export
setMethod("rownames<-", "MatrixList", function(x, value){
  
  rownames(x@unlistData) <- value
  x
})

#' @rdname MatrixList-class
#' @export
setMethod("colnames", "MatrixList", function(x){
  
  colnames(x@unlistData)
  
})

#' @rdname MatrixList-class
#' @export
setMethod("colnames<-", "MatrixList", function(x, value){
  
  colnames(x@unlistData) <- value
  x
})


#' @rdname MatrixList-class
#' @export
setMethod("length", "MatrixList", function(x){
  
  length(x@partitioning)
  
})


#' @rdname MatrixList-class
#' @export
setMethod("width", "MatrixList", function(x){
  
  sapply(x@partitioning, length)
  
})


#' @rdname MatrixList-class
#' @export
setMethod("dim", "MatrixList", function(x){
  
  dim(x@unlistData)
  
})


#' @rdname MatrixList-class
#' @export
setMethod("nrow", "MatrixList", function(x){
  
  nrow(x@unlistData)
  
})

#' @rdname MatrixList-class
#' @export
setMethod("ncol", "MatrixList", function(x){
  
  ncol(x@unlistData)
  
})


#' @rdname MatrixList-class
#' @export
setMethod("[[", "MatrixList", function(x, i){
  
  x@unlistData[x@partitioning[[i]], , drop = FALSE]
  
})

#' @rdname MatrixList-class
#' @param name Name of submatrix to be accessed.
#' @export
setMethod("$", "MatrixList", function(x, name){
  
  x[[name]]
  
})


##############################################################

#' @rdname MatrixList-class
#' @param i,j Indexes of MatrixList to be subsampled.
#' @export
setMethod("[", "MatrixList", function(x, i, j){
  
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







