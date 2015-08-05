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
