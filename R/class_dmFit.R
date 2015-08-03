setClass("dmFit", 
  representation(
    proportions = "MatrixList",
    statistics = "DataFrame"))


setMethod("show", "dmFit", function(object){
  
  cat("\nSlot \"proportions_full\":\n")
  print(object@proportions)
  
  cat("\nSlot \"statistics_full\":\n")
  print(object@statistics)
  
  })

















































