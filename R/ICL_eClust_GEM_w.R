
#' @import igraph
#' @import dplyr
#' @export




evalConditionalLik_wn2 <- function(Z, S1, R1, U, V, Y, 
                                   Eta, W_adj, Lambda_a, phi, EE, 
                                   log_t_k, t_0, alpha,
                                   alpha_0 = 1.1, beta_0 = 25){
  
  
  
  K = nrow(Y)
  M = nrow(EE)
  n = length(S1)
  p = ncol(Y)
  
  
  
  eUiWk = exp(S1 %*% t(rep(1,K)) + U %*% t(Y))
  eViWk = exp(R1 %*% t(rep(1,K)) + V %*% t(Y))
  
  
  fuk = colSums(eUiWk)
  fvk = colSums(eViWk)
  
  A_eta = Eta^2/2
  Pr_w_eta = matrix(1, ncol = K, nrow = M)
  a_phi = phi^2
  W_adj_exp = W_adj %*% t(rep(1,K))
  log_h_phi = -W_adj_exp^2 %*% diag(1/(phi^2))/2 - matrix(1, M, K)%*% diag(log(phi))
  
  
  ret = 0
  
  
  for (m in 1:M){
    if (Z[m] == 1){
      
      ret = ret + log(t_0) - log(n * (n-1)/Lambda_a) -Lambda_a*W_adj[m]
    
      }else{
      
     ret = ret + (log_h_phi[m, (Z[m]-1)] + (Eta[m, (Z[m]-1)] * W_adj[m] - A_eta[m,(Z[m]-1)])/a_phi[(Z[m]-1)] - 
                     
                     log(Pr_w_eta[m,(Z[m]-1)])  +
                     
                     S1[EE[m,1]] + R1[EE[m,2]] + U[EE[m,1],] %*% Y[(Z[m]-1),] + 
                     
                     V[EE[m,2],] %*% Y[(Z[m]-1),] - 
                     
                     log(fuk[(Z[m]-1)]) - log(fvk[(Z[m]-1)] - exp(R1[EE[m,1]] + V[EE[m,1],] %*% Y[(Z[m]-1),])) +
                     
                     log_t_k[Z[m] - 1] + log(1 - t_0))
      
    }
    
    
    
    ret = ret + lgamma(K*alpha) - K*lgamma(alpha) + (alpha_0 - 1)*log(t_0) +
      (beta_0 - 1)*log(1 - t_0) + alpha*sum(log_t_k)
    
  }
  
  return(ret)
  
}

#' @export

ICL.eClust_GEM_w = function(eCl,A){
 
  n = vcount(A)
  
  K = nrow(eCl$estimates$Y)
  
  p = ncol(eCl$estimates$U)
  
  M = ecount(A)
 
  zHat = eCl$estimates$Z_est
  
  n_k = sapply(1:(K+1),function(k) sum(zHat == k))
  
  n_0 = n_k[1]
  
  n_k = n_k[-1]
  
  zHat =as.numeric(as.factor(zHat))
  
  YY = eCl$estimates$Y[n_k>0,]
  
  Elog_t_k <- eCl$estimates$Elog_t_k[n_k>0] 
  
  Eta <- eCl$estimates$Eta[,n_k>0]
  
  phi <- eCl$estimates$phi[n_k>0]
  
  n_k = n_k[n_k>0]
  
  
  ICL = - evalConditionalLik_wn2(zHat,
                                 eCl$estimates$S1,eCl$estimates$R1,  
                                 eCl$estimates$U,eCl$estimates$V,
                                 YY, Eta,
                                 log(E(A)$weight), eCl$userInputs$lambda_a,
                                 phi,
                                 as_edgelist(A,FALSE),
                                 Elog_t_k, 
                                 eCl$estimates$E_t_0, 
                                 eCl$estimates$alpha,
                                 alpha_0 = 1.1, beta_0 = 25)  -
    0.5*log(ecount(A))*( 2^(is_directed(A))*(2*vcount(A)+length(eCl$estimates$U)) + length(YY) +
                           length(phi)*p*p + 2*length(phi) +
                           3*length(eCl$estimates$sigma_SR1) -1 )
  
  
  
  return(c(ICL = ICL))
  
}


