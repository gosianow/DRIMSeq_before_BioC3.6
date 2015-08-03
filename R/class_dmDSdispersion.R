setClassUnion("numericORNULL", c("numeric", "NULL"))

setClass("dmDSdispersion", 
  contains = "dmDSdata",
  representation(mean_expression = "numericORNULL", 
    common_dispersion = "numericORNULL",
    tagwise_dispersion = "numericORNULL"), 
  prototype(mean_expression = NULL,
    common_dispersion = NULL, 
    tagwise_dispersion = NULL))



show_numeric <- function(object, nhead = 2, ntail = 2){
  
  if(!is.null(object)){
    
    nl <- length(object)  
    cat(class(object), "of length", length(object), "\n")
    
    if(nl < (nhead + ntail + 1L)) {
      out <- round(object, 2)
      } else {
        dots <- "... ..."
        names(dots) <- "... ..."
        out <- c(round(head(object, nhead), 2), dots , round(tail(object, ntail), 2))
      }
      print(out, quote = FALSE, right = TRUE)
      
      }else{
        
        print(object)
        
      }

    }




setMethod("show", "dmDSdispersion", function(object){
      
  # cat("Slot \"counts\":\n")
  # print(object@counts)
  
  # cat("Slot \"samples\":\n")
  # print(object@samples)
  
  callNextMethod(object)
  
  cat("\nSlot \"mean_expression\":\n")
  show_numeric(object@mean_expression)
  
  cat("\nSlot \"common_dispersion\":\n")
  show_numeric(object@common_dispersion)
  
  cat("\nSlot \"tagwise_dispersion\":\n")
  show_numeric(object@tagwise_dispersion)
  
  })




setGeneric("dmDSdispersion", function(x, ...) standardGeneric("dmDSdispersion"))

setMethod("dmDSdispersion", "dmDSdata", function(x, mean_expression = TRUE, common_dispersion = FALSE, tagwise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = c("optimize", "optim", "constrOptim", "grid")[4], disp_interval = c(0, 1e+5), disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = c("none", "common", "trended")[1], disp_prior_df = 10, disp_span = 0.3, prop_mode = c( "constrOptim", "constrOptimG", "FisherScoring")[2], prop_tol = 1e-12, verbose = FALSE, BPPARAM = MulticoreParam(workers=1)){
   
  if(mean_expression){
    mean_expression <- dm_estimateMeanExpression(counts = x@counts, BPPARAM = BPPARAM)
    }else{
      mean_expression <- NULL
    }
  
  if(common_dispersion){
    common_dispersion <- dmDS_estimateCommonDispersion(counts = x@counts, samples = x@samples, disp_adjust = disp_adjust, disp_interval = disp_interval, disp_tol = 1e+01, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
    }else{
      common_dispersion <- NULL
    }

  
  if(tagwise_dispersion){
    
    if(disp_mode == "grid" && !is.null(common_dispersion)){
      cat("!Using common dispersion as disp_init in 'grid' mode!\n")
      disp_init <- common_dispersion
    }

    tagwise_dispersion <- dmDS_estimateTagwiseDispersion(counts = x@counts, samples = x@samples, disp_adjust = disp_adjust, disp_mode = disp_mode, disp_interval = disp_interval, disp_tol = disp_tol, disp_init = disp_init, disp_init_weirMoM = disp_init_weirMoM, disp_grid_length = disp_grid_length, disp_grid_range = disp_grid_range, disp_moderation = disp_moderation, disp_prior_df = disp_prior_df, disp_span = disp_span, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
    
  }else{
    tagwise_dispersion <- NULL
  }
  
  
  return(new("dmDSdispersion", mean_expression = mean_expression, common_dispersion = common_dispersion, tagwise_dispersion = tagwise_dispersion, counts = x@counts, samples = x@samples))
  
  
  })



setMethod("dmDSdispersion", "dmDSdispersion", function(x, mean_expression = FALSE, common_dispersion = FALSE, tagwise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = c("optimize", "optim", "constrOptim", "grid")[4], disp_interval = c(0, 1e+5), disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = c("none", "common", "trended")[1], disp_prior_df = 10, disp_span = 0.3, prop_mode = c( "constrOptim", "constrOptimG", "FisherScoring")[2], prop_tol = 1e-12, verbose = FALSE, BPPARAM = MulticoreParam(workers=1)){
   
  if(mean_expression){
    mean_expression <- dm_estimateMeanExpression(counts = x@counts, BPPARAM = BPPARAM)
    }else{
      mean_expression <- x@mean_expression
    }
  
  if(common_dispersion){
    common_dispersion <- dmDS_estimateCommonDispersion(counts = x@counts, samples = x@samples, disp_adjust = disp_adjust, disp_interval = disp_interval, disp_tol = 1e+01, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
    }else{
      common_dispersion <- x@common_dispersion
    }

  
  if(tagwise_dispersion){
    
    if(disp_mode == "grid" && !is.null(common_dispersion)){
      cat("!Using common_dispersion =", round(common_dispersion, 2), "as disp_init in 'grid' mode!\n")
      disp_init <- common_dispersion
    }

    tagwise_dispersion <- dmDS_estimateTagwiseDispersion(counts = x@counts, samples = x@samples, disp_adjust = disp_adjust, disp_mode = disp_mode, disp_interval = disp_interval, disp_tol = disp_tol, disp_init = disp_init, disp_init_weirMoM = disp_init_weirMoM, disp_grid_length = disp_grid_length, disp_grid_range = disp_grid_range, disp_moderation = disp_moderation, disp_prior_df = disp_prior_df, disp_span = disp_span, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
    
  }else{
    tagwise_dispersion <- x@tagwise_dispersion
  }
  
  
  return(new("dmDSdispersion", mean_expression = mean_expression, common_dispersion = common_dispersion, tagwise_dispersion = tagwise_dispersion, counts = x@counts, samples = x@samples))
  
  
  })




setGeneric("dmDSplotDispersion", function(x, ...) standardGeneric("dmDSplotDispersion"))

setMethod("dmDSplotDispersion", "dmDSdispersion", function(x, out_dir = NULL){
  
  # tagwise_dispersion = x@tagwise_dispersion; mean_expression = x@mean_expression; nr_features = width(x@counts@partitioning); common_dispersion = x@common_dispersion
  
  dmDS_plotDispersion(tagwise_dispersion = x@tagwise_dispersion, mean_expression = x@mean_expression, nr_features = width(x@counts@partitioning), common_dispersion = x@common_dispersion, out_dir = out_dir)
  
  })












































