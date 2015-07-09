show.dmDSdata <- function(object){
  
  cat("@counts\n", sep="")
  
  x <- object@counts
  n <- length(x)
  if(n) {
    i <- names(x)
    if(is.null(i)) i <- seq_len(n)
    if(length(i) >= 2) i <- i[1:2]
    
    for (what in i) {
      y <- x[[what]]
      cat("$", what, "\n", sep="")
      
      d <- dim(y)
      if(!is.null(d)){
        
        if(any(d > 10)) {
          
          if(d[1] > 10)
            y <- y[1:6,]
          
          if(d[2] > 10)
            y <- y[, 1:6]
          
          print(y)
          
          if(d[1] > 10)
            cat("*",d[1] - 6,"more rows ...\n") 
          
          if(d[2] > 10)
            cat("*",d[2] - 6, "more columns ...\n")
          
        } else
          print(y)
        
      } else
        print(y) 
      
    }
    cat("**", n-2, "more elements ...\n\n")
  }
  
  cat("@samples\n", sep="")
  x <- object@samples
  d <- dim(x)
  n <- d[1]
  if(n > 10) {
    
    x <- x[1:6, , drop = FALSE]
    
    if(d[2] > 10)
      x <- x[, 1:6, drop = FALSE]
    
    print(x)
    cat("*",n-6,"more rows ...\n") 
    
    if(d[2] > 10)
      cat("*",d[2] - 6, "more columns ...\n")
    
  } else
    print(x)
  
}

setMethod("show", signature(object = "dmDSdata"), show.dmDSdata)



show.dmSQTLdata <- function(object){
  
  cat("@counts\n", sep="")
  
  x <- object@counts
  n <- length(x)
  if(n) {
    i <- names(x)
    if(is.null(i)) i <- seq_len(n)
    if(length(i) >= 2) i <- i[1:2]
    
    for (what in i) {
      y <- x[[what]]
      cat("$", what, "\n", sep="")
      
      d <- dim(y)
      if(!is.null(d)){
        
        if(any(d > 10)) {
          
          if(d[1] > 10)
            y <- y[1:6,]
          
          if(d[2] > 10)
            y <- y[, 1:6]
          
          print(y)
          
          if(d[1] > 10)
            cat("*",d[1] - 6,"more rows ...\n") 
          
          if(d[2] > 10)
            cat("*",d[2] - 6, "more columns ...\n")
          
        } else
          print(y)
        
      } else
        print(y) 
      
    }
    cat("**", n-2, "more elements ...\n\n")
  }
  
  
  cat("@genotypes\n", sep="")
  
  x <- object@genotypes
  n <- length(x)
  if(n) {
    i <- names(x)
    if(is.null(i)) i <- seq_len(n)
    if(length(i) >= 2) i <- i[1:2]
    
    for (what in i) {
      y <- x[[what]]
      cat("$", what, "\n", sep="")
      
      d <- dim(y)
      if(!is.null(d)){
        
        if(any(d > 10)) {
          
          if(d[1] > 10)
            y <- y[1:6,]
          
          if(d[2] > 10)
            y <- y[, 1:6]
          
          print(y)
          
          if(d[1] > 10)
            cat("*",d[1] - 6,"more rows ...\n") 
          
          if(d[2] > 10)
            cat("*",d[2] - 6, "more columns ...\n")
          
        } else
          print(y)
        
      } else
        print(y) 
      
    }
    cat("**", n-2, "more elements ...\n\n")
  }
  
  
  cat("@samples\n", sep="")
  x <- object@samples
  d <- dim(x)
  n <- d[1]
  if(n > 10) {
    
    x <- x[1:6, , drop = FALSE]
    
    if(d[2] > 10)
      x <- x[, 1:6, drop = FALSE]
    
    print(x)
    cat("*",n-6,"more rows ...\n") 
    
    if(d[2] > 10)
      cat("*",d[2] - 6, "more columns ...\n")
    
  } else
    print(x)
  
  
}

setMethod("show", signature(object = "dmSQTLdata"), show.dmSQTLdata)

















