
dm_cpm <- function(x){

  lib_size <- 1e-6 * colSums(x, na.rm = TRUE)
  t(t(x)/lib_size)
  
}

