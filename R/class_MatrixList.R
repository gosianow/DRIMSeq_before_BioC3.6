setClassUnion("matrixORNULL", c("matrix", "NULL"))

#' @importClassesFrom IRanges CompressedList
setClass("MatrixList", contains = "CompressedList", representation(unlistData = "matrix"), prototype(elementType = "matrix"))

### Constructor:

MatrixList <- function(...){
  
  listData <- list(...)
  if (length(listData) == 1L && is.list(listData[[1L]]))
    listData <- listData[[1L]]
  if (length(listData) == 0L) {
    return(new("MatrixList"))
  } else {
    if (!all(sapply(listData, is, "matrixORNULL")))
      stop("all elements in '...' must be matrix or NULL objects")
    
    listData <- listData[!sapply(listData, is.null)]
    
    unlistData <- do.call(rbind, listData)
    
    return(new("MatrixList", unlistData = unlistData, partitioning = IRanges::PartitioningByEnd(listData)))
    
  }
}

# nhead = 3; ntail = 3

show_matrix <- function(object, nhead = 3, ntail = 3){
  nr <- nrow(object)
  nc <- ncol(object)
  
  cat(mode(object), " " ,class(object), " with ", nr, ifelse(nr == 1, " row and ", " rows and "), nc, ifelse(nc == 1, " column\n", " columns\n"), sep = "")
  
  if(nr > 0 && nc > 0){
    nms <- rownames(object)
    if(nr < (nhead + ntail + 1L)) {
      out <- object
    } else {
      
      out <- do.call(rbind, list(head(object, nhead), rep.int("...", nc), tail(object, ntail)))
      
      rownames(out) <- rownames_matrix(nms, nr, nhead, ntail)
    }
    
    ### print adjusted for numeric
    print(out, quote = FALSE, right = TRUE, na.print = "NA")
  }
  
}


rownames_matrix <- function(nms, nrow, nhead, ntail){
  p1 <- ifelse (nhead == 0, 0L, 1L)
  p2 <- ifelse (ntail == 0, 0L, ntail-1L)
  s1 <- s2 <- character(0)
  
  if (is.null(nms)) {
    if (nhead > 0)
      s1 <- paste0("[",as.character(p1:nhead), ",]")
    if (ntail > 0)
      s2 <- paste0("[",as.character((nrow-p2):nrow), ",]")
  } else {
    if (nhead > 0)
      s1 <- paste0(head(nms, nhead))
    if (ntail > 0)
      s2 <- paste0(tail(nms, ntail))
  }
  c(s1, "...", s2)
}


setMethod("show", "MatrixList", function(object){
  
  nn <- length(object)
  cat("MatrixList of length", nn,"\n")
  
  if(nn){
    
    i <- names(object)
    if(is.null(i)) i <- seq_len(nn)
    
    what <- i[1]
    cat("$", what, "\n", sep="")
    show_matrix(object@unlistData[object@partitioning[[what]], , drop = FALSE])
    
    if(nn > 1){
      if(nn > 2)
        cat(":\n", sep="")
      what <- i[nn]
      cat("$", what, "\n", sep="")
      show_matrix(object@unlistData[object@partitioning[[what]], , drop = FALSE])
    }
    
  }
  
  if(!is.null(object@elementMetadata))
  print(object@elementMetadata)
  
})



setMethod("[[", "MatrixList", function(x, i, ...){
  
  # di <- NULL 
  # try(di <- x@unlistData[x@partitioning[[i]], , drop = FALSE], silent = TRUE) 
  # di
  
  x@unlistData[x@partitioning[[i]], , drop = FALSE]
  
})




setMethod("c", "MatrixList", function(x, ..., recursive = FALSE){
  
  if (!identical(recursive, FALSE))
    stop("\"c\" method for MatrixList objects ", "does not support the 'recursive' argument")
  
  if (missing(x))
    tls <- unname(list(...))
  else tls <- unname(list(x, ...))
  if (!all(sapply(tls, is, "MatrixList")))
    stop("all arguments in '...' must be MatrixList objects")
  ecs <- sapply(tls, elementType)
  if (!all(sapply(ecs, extends, ecs[[1L]])))
    stop("all arguments in '...' must have an element class that extends that of the first argument")
  
  
  unlistData <- do.call(rbind, lapply(tls, function(tl) tl@unlistData))
  
  ends <- cumsum(do.call(c, lapply(tls, function(tl){
    # tl <- tls[[1]]
    
    w <- IRanges::width(tl@partitioning)
    if(!is.null(names(tl@partitioning)))
      names(w) <- names(tl@partitioning)
    w
    
  })))
  
  if(all(sapply(tls, function(tl) is.null(tl@elementMetadata))))
    elementMetadata <- NULL
  
  dim_meta <- unique(unlist(lapply(tls, function(tl) ncol(tl@elementMetadata))))
  
  if(!length(dim_meta) == 1){
    elementMetadata <- NULL
  }else{
    
    names_meta <- do.call(rbind, lapply(tls, function(tl) colnames(tl@elementMetadata)))
    
    if(!nrow(unique(names_meta)) == 1){
      elementMetadata <- NULL
      
    }else{
      elementMetadata <- do.call(rbind, lapply(tls, function(tl){
        # tl <- tls[[3]]
        if(is.null(tl@elementMetadata)){
          el_meta <- S4Vectors::DataFrame(matrix(NA, length(tl@partitioning), dim_meta))
          colnames(el_meta) <- names_meta[1, ]
          # if(!is.null(names(tl@partitioning)))
          #   rownames(el_meta) <- names(tl@partitioning)
          return(el_meta)
        }
        tl@elementMetadata
      }))
    }
  }
  
  new("MatrixList", unlistData = unlistData, partitioning = IRanges::PartitioningByEnd(ends), elementMetadata = elementMetadata)
  
})



setMethod("nrow", "MatrixList", function(x){
  
  nrow(x@unlistData)
  
})


setMethod("ncol", "MatrixList", function(x){
  
  ncol(x@unlistData)
  
})






