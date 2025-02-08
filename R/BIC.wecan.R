#' BIC for wecan objects
#' 
#' @param object object of class 'wecan'
#' 
#' @import igraph
#' @import dplyr
#' @export BIC.wecan
#' @exportS3Method WECAN::BIC


BIC.wecan = function(object){
  c(BIC = -2*evalMargLogLik_w(object$estimates$S1,object$estimates$R1,  
                              
                              object$estimates$U,object$estimates$V,
                              
                              object$estimates$Y, object$estimates$Eta,
                              
                              E(object$A)$weight, object$estimates$alpha,
                              
                              object$estimates$phi,
                              
                              as_edgelist(object$A,FALSE),
                              
                              object$estimates$Pmk) +
      
      log(ecount(object$A))*( 2^(is_directed(object$A))*(2*vcount(object$A)+length(object$estimates$U)) + # S1,S2,R1,R2,U,V
                         
                         length(object$estimates$Y) + # Y
                         
                         length(object$estimates$Lambda)+ length(object$estimates$phi) + # Lambda, phi
                         
                         2*length(object$estimates$sigma_SR1) +  length(object$estimates$sigma_UV) + #sigma_SR1, sigma_SR2, sigma_UV
                         
                         2*length(object$estimates$alpha) -1 ))  # alpha, beta
  
}

