##############################################################################
# calculate tagwise dispersions 
# dmEstimateTagwiseDisp, dmSQTLEstimateTagwiseDisp
##############################################################################

## Computes the MoM estimate of theta (and std. error) / from dirmult package / use as starting values for gamma0 = (1-mom)/mom

weirMoM <- function(data, se=FALSE){

#   data <- t(data)
  
  K <- ncol(data)
  J <- nrow(data)
  MoM <- colSums(data)/sum(data)
  Sn <- rowSums(data)
  MSP <- (J-1)^(-1)*sum(rowSums((data/rowSums(data)-matrix(rep(MoM,J),J,K,byrow=T))^2)*Sn)
  MSG <- (sum(data)-J)^(-1)*sum(rowSums(data/rowSums(data)*(1-data/rowSums(data)))*Sn)
  nc <- 1/(J-1)*(sum(Sn)-sum(Sn^2)/sum(Sn))
 
  MoM.wh <- (MSP-MSG)/(MSP+(nc-1)*MSG)
  
#   if(se){
#     ## Formula by Li, ref in Weir-Hill 2002
#     std.er <- sqrt(2*(1-MoM.wh)^2/(J-1)*((1+(nc-1)*MoM.wh)/nc)^2)
#     list(theta=MoM.wh,se=std.er)
#   }
#   else MoM.wh
  
  if(MoM.wh <= 0)
    MoM.wh <- 0.1 # MoM.wh <- 0.005
  
  return((1-MoM.wh)/MoM.wh)
  
}

