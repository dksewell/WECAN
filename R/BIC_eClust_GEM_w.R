
#' @import igraph
#' @import dplyr
#' @export



BIC_eClust_GEM_w = function(eCl,A){
  c(BIC = -2*evalMargLogLik_w(eCl$estimates$S1,eCl$estimates$R1,  
                              
                              eCl$estimates$U,eCl$estimates$V,
                              
                              eCl$estimates$Y, eCl$estimates$Eta,
                              
                              log(E(A)$weight), eCl$estimates$alpha,
                              
                              eCl$estimates$phi,
                              
                              as_edgelist(A,FALSE),
                              
                              eCl$estimates$Pmk) +
      
      log(ecount(A))*( 2^(is_directed(A))*(2*vcount(A)+length(eCl$estimates$U)) + # S1,S2,R1,R2,U,V
                         
                         length(eCl$estimates$Y) + # Y
                         
                         length(eCl$estimates$Lambda)+ length(eCl$estimates$phi) + # Lambda, phi
                         
                         2*length(eCl$estimates$sigma_SR1) +  length(eCl$estimates$sigma_UV) + #sigma_SR1, sigma_SR2, sigma_UV
                         
                         2*length(eCl$estimates$alpha) -1 ))  # alpha, beta
  
}

