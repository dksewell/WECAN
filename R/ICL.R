#' Integrated Complete Likelihood (ICL) for wecan objects
#' 
#' @param object Object of class wecan
#' 
#' @import igraph
#' @import dplyr
#' @export ICL

ICL = function(object){
 
  n = vcount(object$A)
  
  K = nrow(object$estimates$Y)
  
  p = ncol(object$estimates$U)
  
  M = ecount(object$A)
 
  zHat = object$estimates$Z_est
  
  n_k = sapply(1:(K+1),function(k) sum(zHat == k))
  
  n_0 = n_k[1]
  
  n_k = n_k[-1]
  
  zHat =as.numeric(as.factor(zHat))
  
  YY = object$estimates$Y[n_k>0,]
  
  Elog_t_k <- object$estimates$Elog_t_k[n_k>0] 
  
  Eta <- object$estimates$Eta[,n_k>0]
  
  phi <- object$estimates$phi[n_k>0]
  
  n_k = n_k[n_k>0]
  
  
  ICL = - evalConditionalLik_wn2(zHat,
                                 object$estimates$S1,
                                 object$estimates$R1,  
                                 object$estimates$U,
                                 object$estimates$V,
                                 YY, 
                                 Eta,
                                 igraph::edge_attr(object$A,object$userInputs$weight_attr_name),
                                 object$userInputs$lambda_a,
                                 phi,
                                 as_edgelist(object$A,FALSE),
                                 Elog_t_k, 
                                 object$estimates$E_t_0, 
                                 object$estimates$alpha,
                                 alpha_0 = 1.1, beta_0 = 25)  -
    0.5*log(ecount(object$A))*( 2^(is_directed(object$A))*(2*vcount(object$A)+length(object$estimates$U)) + length(YY) +
                           length(phi)*p*p + 2*length(phi) +
                           3*length(object$estimates$sigma_SR1) -1 )
  
  
  
  return(c(ICL = ICL))
  
}


