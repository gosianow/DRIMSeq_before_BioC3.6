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
            y <- y[1:5,]
          
          if(d[2] > 10)
            y <- y[, 1:5]
          
          print(y)
          
          if(d[1] > 10)
            cat("*",d[1] - 5,"more rows ...\n") 
          
          if(d[2] > 10)
            cat("*",d[2] - 5, "more columns ...\n")
          
        } else
          print(y)
        
      } else
        print(y) 
      
    }
    cat("**", n-2, "more elements ...\n\n")
  }
	
  cat("@samples\n", sep="")
  print(object@samples)
  
}

setMethod("show", signature(object = "dmDSdata"), show.dmDSdata)