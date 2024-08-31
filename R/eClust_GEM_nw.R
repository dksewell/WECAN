
#' @import igraph
#' @import MASS
#' @import LaplacesDemon
#' @import RcppArmadillo
#' @import Rcpp
#' @import aricode
#' @import aLSEC
#' @import mclust
#' @import aLSEC
#' @export



eClust_GEM_wn = function(A, K = 2,
                         p = 3,
                         maxIter=1e3,
                         maxIterVB=100,
                         a_0 = 1,
                         A_0 = 3,
                         B_0 = 0.5,
                         nu_0 = 2,
                         eta_0 = 5,
                         QNSteps=25,
                         CGSteps=25,
                         distribution = "Normal",
                         if_log = FALSE,
                         eps=1e-5,
                         upsilon_SR1 = 2,
                         upsilon_SR2 = 2,
                         upsilon_UV = 2,
                         Phi_SR1 = toeplitz((2:1)/2),
                         Phi_SR2 = toeplitz((2:1)/2),
                         Phi_UV = toeplitz((2:1)/2),
                         noise_mean = 0.1,
                         alpha_init = NULL,
                         a_a = 1,
                         b_a = 200,
                         t_0 = 0.05,
                         alpha_0 = 1.1,
                         beta_0 = 25){
  
  
  ### Create objects
  M = ecount(A)
  n = vcount(A)
  EE = cbind(ends(A,1:ecount(A),FALSE))
  W = E(A)$weight
  A_adj <- as_adjacency_matrix(A, attr = "weight")
  
  
  ### sourceCpp(file = paste0(wd, "LSEC/src/indexEdges.cpp"))
  cat("\n Indexing edgelist \n")
  temp = indexEdges(EE,n)
  Mi1 = temp$Mi1[,1:max(temp$Mi1Index)]
  Mi2 = temp$Mi2[,1:max(temp$Mi2Index)]
  Mi1Index = temp$Mi1Index
  Mi2Index = temp$Mi2Index
  
  rm(temp)
  
  ## --- Helper functions:
  test_fun = function(x){
    QUVY(SS1 = x[1:n],
         SS2 = x[(n+1):(2*n)],
         U = matrix(x[(2*n + 1):(2*n + n*p)],n,p),
         RR1 = x[(2*n + n*p + 1):(3*n + n*p)],
         RR2 = x[(3*n + n*p + 1):(4*n + n*p)],
         V = matrix(x[(4*n + n*p + 1):(4*n + 2*n*p)],n,p),
         Y = matrix(x[(4*n + 2*n*p + 1):(4*n + 2*n*p + K*p)],K,p),
         Lambda =  array(x[(4*n + 2*n*p + K*p + 1):(4*n + 2*n*p + K*p + K*p*p)],
                         dim = c(p,p,K)),
         beta = x[(4*n + 2*n*p + K*p + K*p*p + 1):(4*n + 2*n*p + K*p + K*p*p + K)],
         phi = x[(length(x) - K + 1):length(x)],
         lam = lam, W_adj, Pmk, Pk,EE,
         Eta = computeEta(SS2 = x[(n+1):(2*n)], RR2 = x[(3*n + n*p + 1):(4*n + n*p)],
                          beta = x[(4*n + 2*n*p + K*p + K*p*p + 1):(4*n + 2*n*p + K*p + K*p*p + K)],
                          U = matrix(x[(2*n + 1):(2*n + n*p)],n,p),
                          V = matrix(x[(4*n + n*p + 1):(4*n + 2*n*p)],n,p),
                          EE,
                          Lambda = array(x[(4*n + 2*n*p + K*p + 1):(4*n + 2*n*p + K*p + K*p*p)],
                                         dim = c(p,p,K)), K),
         sigma_SR1 = sigma_SR1,
         sigma_SR2 = sigma_SR2,sigma_UV = sigma_UV,
         SR1 = cbind(x[1:n],  x[(2*n + n*p + 1):(3*n + n*p)]),
         SR2 = cbind(x[(n+1):(2*n)], x[(3*n + n*p + 1):(4*n + n*p)]),
         UV = cbind(matrix(x[(2*n + 1):(2*n + n*p)],n,p),
                    matrix(x[(4*n + n*p + 1):(4*n + 2*n*p)],n,p)),
         nu_0 = nu_0, eta_0 = eta_0)
  }
  
  test_gr = function(x){
    dQUVY(SS1 = x[1:n],
          SS2 = x[(n+1):(2*n)],
          U = matrix(x[(2*n + 1):(2*n + n*p)],n,p),
          RR1 = x[(2*n + n*p + 1):(3*n + n*p)],
          RR2 = x[(3*n + n*p + 1):(4*n + n*p)],
          V = matrix(x[(4*n + n*p + 1):(4*n + 2*n*p)],n,p),
          Y = matrix(x[(4*n + 2*n*p + 1):(4*n + 2*n*p + K*p)],K,p),
          Lambda = array(x[(4*n + 2*n*p + K*p + 1):(4*n + 2*n*p + K*p + K*p*p)],
                         dim = c(p,p,K)),
          beta = x[(4*n + 2*n*p + K*p + K*p*p + 1):(4*n + 2*n*p + K*p + K*p*p + K)],
          phi = x[(length(x) - K + 1):length(x)],
          lam, W_adj, Pmk, Pmki1, Pmki2,
          Pk,Mi1Index,Mi2Index,Mi1,Mi2,
          EE,  
          Eta = computeEta(SS2 = x[(n+1):(2*n)], RR2 = x[(3*n + n*p + 1):(4*n + n*p)],
                           beta = x[(4*n + 2*n*p + K*p + K*p*p + 1):(4*n + 2*n*p + K*p + K*p*p + K)],
                           U = matrix(x[(2*n + 1):(2*n + n*p)],n,p),
                           V = matrix(x[(4*n + n*p + 1):(4*n + 2*n*p)],n,p),
                           EE,
                           Lambda = array(x[(4*n + 2*n*p + K*p + 1):(4*n + 2*n*p + K*p + K*p*p)],
                                          dim = c(p,p,K)), K),
          solve(sigma_SR1),solve(sigma_SR2),
          solve(sigma_UV), nu_0 = nu_0, eta_0 = eta_0)
  }
  
  
  fun_helper <- function(x){x %*% t(x)}
  
  find_clust <- function(Pmki){
    return(order(Pmki, decreasing = TRUE)[1])
  }
  
  
  #--- End helper functions
  
  ### Initialize
  
  Lambda= array(0,dim = c(p,p,K))
  
  beta = rep(0,K)
  
  beta_int = beta
  
  
  # S1 = degree(A, mode = "out")
  # S1 = (S1 - mean(S1))/max(abs(S1 - mean(S1)))
  # S1_init = S1
  # R1 = degree(A, mode = "in")
  # R1 = (R1 - mean(R1))/max(abs(R1 - mean(R1)))
  # R1_init = R1
  # SR1 = cbind(S1,R1)
  
  
  S2 = strength(A, mode = "out")
  S2 = (S2 - mean(S2))/max(abs(S2 - mean(S2)))
  S2_init = S2
  R2 = strength(A, mode = "in")
  R2 = (R2 - mean(R2))/max(abs(R2 - mean(R2)))
  R2_init = R2
  SR2 = cbind(S2,R2)
  
  
  # U = svd(A_adj,nu = p, nv = p)$u
  # V = svd(A_adj,nu = p, nv = p)$v
  # 
  # U_init = U
  # V_init = V
  # UV = cbind(U,V)
  # 
  # Y = matrix(1,K,p)
  # alpha = ifelse(is.null(alpha_init), a_a/b_a, alpha_init)
  
  lam = 0.1
  
  sigma_SR1 = rinvwishart(upsilon_SR1, Phi_SR1*0.1)
  sigma_SR2 = rinvwishart(upsilon_SR2, Phi_SR2*0.1)
  sigma_UV = rinvwishart(upsilon_UV, Phi_UV*0.1)
  lambda_a = 1/noise_mean
  
  
  
  #############################################################################
  cat("\n Finding reasonable starting S1 R1 U V and Y \n")
  
  noise_index <- which(W <= quantile(W, t_0))
  non_noise_index <- setdiff(1:M, noise_index)
  A2 <- delete.edges(A, noise_index)
  
  init <- eClustaLSEC_GEM(A2,K=K,p=p)
  
  
  S1 <- init$estimates$S
  R1 <- init$estimates$R
  SR1 = cbind(S1,R1)
  
  
  U <- init$estimates$U
  V <- init$estimates$V
  UV = cbind(U,V)
  
  
  Y <- init$estimates$W
  alpha <- init$estimates$alpha
  
  # plot(U, vegan::procrustes(Truth$U, U)$Yrot)
  # plot(V, vegan::procrustes(Truth$V, V)$Yrot)
  
  
  
  
  Eta = computeEta(S2, R2, beta, U, V, EE, Lambda, K)
  
  if (distribution == "Poisson"){
    W_adj = ceiling(W)
    h_w_phi = 1/factorial(W_adj)
    A_eta = exp(Eta) 
    Pr_w_eta = 1 - exp(-exp(Eta))
    A_eta_dev = exp(Eta) 
    Pr_eta_dev = exp(-exp(Eta))*exp(Eta)
    
    phi = rep(0,K)
    get_pi_phi <- function(phi) {1}
    a_phi = rep(1,K)
    a_phi_dev = rep(1,K)
    pi_phi = get_pi_phi(phi)
  }
  
  if (distribution == "Normal"){
    if (if_log) {W_adj = log(W)
    }else {W_adj = W}
    
    phi = rep(sd(W_adj), K)
    # phi = rep(0.1, K)
    
    A_eta = Eta^2/2
    Pr_w_eta = matrix(1, ncol = K, nrow = M)
    A_eta_dev = Eta
    Pr_eta_dev = matrix(0, ncol = K, nrow = M)
    
    a_phi = phi^2
    a_phi_dev = 2*phi
    
    pi_phi = (1 + phi^2/(nu_0*eta_0^2))^(0.5*(nu_0 + 1))
    
    pi_phi_dev = (nu_0 + 1)*phi/(nu_0*eta_0^2 + phi^2)
    
    
    W_adj_exp = W_adj %*% t(rep(1,K))
    
    log_h_phi = - W_adj_exp^2 %*% diag(1/(phi^2))/2 - matrix(1, M, K)%*% diag(log(phi))
    h_phi = exp(log_h_phi)
    
    h_phi_dev = W_adj_exp^2 %*% diag(1/(phi^3)) - matrix(1, M, K)%*% diag(1/phi)
  }
  
  
  ## sourceCpp(file = "H:/Project/Thesis/Thesis/weighted_edge_cluster/Scripts/Cpp/computePmk_w.cpp")
  Pm0 = numeric(M)
  Pm0[noise_index] = 1
  Pm0[non_noise_index] = 0
  Pmk = init$estimates$Pmk
  
  
  Pmk_all = matrix(0, nrow = M, ncol = (K+1))
  Pmk_all[,1] =  Pm0
  Pmk_all[non_noise_index,2:(K+1)] = Pmk
  
  Pmk = Pmk_all[,2:(K+1)]
  
  
  # alphaTld = rep(alpha, K)
  # Elog_t_k = digamma(alphaTld) - digamma(K*alpha + M - t_0)
  # E_t_0 = (alpha_0+ t_0)/(alpha_0+beta_0 + M)
  # Pmk_all = computePmk_wn(S1,R1,U,V,Y,EE,W_adj,Eta,A_eta,a_phi,Pr_w_eta,h_phi, lambda_a, Elog_t_k, E_t_0)
  # Pmk = Pmk_all[,2:(K+1)]
  # Pm0 = Pmk_all[,1]
  
  
  Pk = colSums(Pmk)
  alphaTld = alpha + colSums(Pmk)
  x_0 = sum(Pm0)
  
  
  Z_prev <- apply(Pmk_all,1,find_clust)
  Z_init <- Z_prev
  
  
  cat("\n Beginning EM algorithm \n")
  pb = txtProgressBar(0,maxIter-1,style=3)
  for(it in 2:maxIter){
    
    ### Update Pmk and Pk (E-step)
    
    Eta = computeEta(S2, R2, beta, U, V, EE, Lambda, K)
    
    if (distribution == "Poisson"){
      W_adj = ceiling(W)
      h_w_phi = 1/factorial(W_adj)
      A_eta = exp(Eta) 
      Pr_w_eta = 1 - exp(-exp(Eta))
      A_eta_dev = exp(Eta) 
      Pr_eta_dev = exp(-exp(Eta))*exp(Eta)
      
      phi = rep(0,K)
      get_pi_phi <- function(phi) {1}
      a_phi = rep(1,K)
      a_phi_dev = rep(1,K)
      pi_phi = get_pi_phi(phi)
    }
    
    if (distribution == "Normal"){
      A_eta = Eta^2/2
      Pr_w_eta = matrix(1, ncol = K, nrow = M)
      A_eta_dev = Eta
      Pr_eta_dev = matrix(0, ncol = K, nrow = M)
      
      a_phi = phi^2
      a_phi_dev = 2*phi
      
      pi_phi = (1 + phi^2/(nu_0*eta_0^2))^(0.5*(nu_0 + 1))
      
      pi_phi_dev = (nu_0 + 1)*phi/(nu_0*eta_0^2 + phi^2)
      
      
      W_adj_exp = W_adj %*% t(rep(1,K))
      
      log_h_phi = - W_adj_exp^2 %*% diag(1/(phi^2))/2 - matrix(1, M, K)%*% diag(log(phi))
      h_phi = exp(log_h_phi)
      
      h_phi_dev = W_adj_exp^2 %*% diag(1/(phi^3)) - matrix(1, M, K)%*% diag(1/phi)
    }
    
    
    Elog_t_k = digamma(alphaTld) - digamma(K*alpha + M - x_0)
    E_t_0 = (alpha_0+ x_0)/(alpha_0+beta_0 + M)
    
    ## sourceCpp(file = "H:/Project/Thesis/Thesis/weighted_edge_cluster/Scripts/Cpp/computePmk_wn.cpp")
    Pmk_all = computePmk_wn(S1,R1,U,V,Y,EE,W_adj,Eta,A_eta,a_phi,Pr_w_eta,h_phi, lambda_a, Elog_t_k, E_t_0)
    Pmk = Pmk_all[,2:(K+1)]
    Pm0 = Pmk_all[,1]
    
    Pk = colSums(Pmk)
    PmklogPmk = sum(Pmk_all[Pmk_all!=0]*log(Pmk_all[Pmk_all!=0]))
    
    ## sourceCpp(file = "H:/Project/Thesis/Thesis/LSEC/src/getPmki.cpp")
    Pmki1 = getPmki(Pmk,Mi1Index,Mi1)
    Pmki2 = getPmki(Pmk,Mi2Index,Mi2)
    
    
    alphaTld = alpha + colSums(Pmk)
    x_0 = sum(Pm0)
    
    
    ELBO_Estep = numeric(maxIterVB)
    
    
    ## sourceCpp(file = "H:/Project/Thesis/Thesis/weighted_edge_cluster/Scripts/Cpp/evalMargPost_w.cpp")
    ELBO_Estep[1] = computeELBO_wn(S1,R1,beta,SR1,SR2,UV,U,V,Y,Eta,log_h_phi,A_eta,a_phi,
                                   pi_phi,Pr_w_eta,  EE,
                                   Phi_SR1,Phi_SR2,Phi_UV,W_adj,a_0,Lambda,lam, A_0, B_0,Pmk, Pk,Pm0,
                                   Elog_t_k, E_t_0, lambda_a, alpha, PmklogPmk, alphaTld, x_0)
    
    
    
    
    for(itVB in 2:maxIterVB){
      # Update alphaTilde; x_0
      alphaTld_old = alphaTld
      alphaTld = alpha + colSums(Pmk)
      x_0_old = x_0
      x_0 = sum(Pm0)
      
      # Update Pmk; Pm0
      Pmk_old = Pmk; Pm0_old = Pm0
      Elog_t_k = digamma(alphaTld) - digamma(K*alpha + M - x_0)
      E_t_0 = (alpha_0+x_0)/(alpha_0 + beta_0 + M)
      
      Pmk_all = computePmk_wn(S1,R1,U,V,Y,EE,W_adj,Eta,A_eta,a_phi,Pr_w_eta,h_phi, lambda_a, Elog_t_k, E_t_0)
      Pmk = Pmk_all[,2:(K+1)]
      Pm0 = Pmk_all[,1]
      PmklogPmk = sum(Pmk_all[Pmk_all!=0]*log(Pmk_all[Pmk_all!=0]))
      
      # Check for convergence
      ELBO_Estep[itVB] = computeELBO_wn(S1,R1,beta,SR1,SR2,UV,U,V,Y,Eta,log_h_phi,A_eta,a_phi,
                                        pi_phi,Pr_w_eta,  EE,
                                        Phi_SR1,Phi_SR2,Phi_UV,W_adj,a_0,Lambda,lam, A_0, B_0,Pmk, Pk,Pm0,
                                        Elog_t_k, E_t_0, lambda_a, alpha, PmklogPmk, alphaTld, x_0)
      
      if( ((ELBO_Estep[itVB] - ELBO_Estep[itVB-1])/abs(ELBO_Estep[itVB-1]) < eps ) &
          ( ELBO_Estep[itVB] <= ELBO_Estep[itVB-1] ) ) break
    }
    
    
    Pk = colSums(Pmk)
    Pmki1 = getPmki(Pmk,Mi1Index,Mi1)
    Pmki2 = getPmki(Pmk,Mi2Index,Mi2)
    
    
    Z_cur = apply(Pmk_all,1,find_clust)
    
    
    ###Check for convergence
    
    if (adjustedRandIndex(Z_cur, Z_prev) >= 0.99) break
    Z_prev = Z_cur
    
    
    ### Update parameters
    Opt = optim(par=c(S1,S2,U,R1,R2,V,Y,Lambda,beta,phi),
                fn = test_fun,
                gr = test_gr,
                method="CG",
                control = list(maxit=100,fnscale=-1,trace=FALSE))
    
    S1 = Opt$par[1:n]
    S2 = Opt$par[(n+1):(2*n)]
    U = matrix(Opt$par[(2*n + 1):(2*n + n*p)],n,p)
    
    R1 = Opt$par[(2*n + n*p + 1):(3*n + n*p)]
    R2 = Opt$par[(3*n + n*p + 1):(4*n + n*p)]
    V = matrix(Opt$par[(4*n + n*p + 1):(4*n + 2*n*p)],n,p)
    
    Y = matrix(Opt$par[(4*n + 2*n*p + 1):(4*n + 2*n*p + K*p)],K,p)
    Lambda = array(Opt$par[(4*n + 2*n*p + K*p + 1):(4*n + 2*n*p + K*p + K*p*p)],
                   dim = c(p,p,K))
    beta = Opt$par[(4*n + 2*n*p + K*p + K*p*p + 1):(4*n + 2*n*p + K*p + K*p*p + K)]
    phi = Opt$par[(length(Opt$par) - K + 1):length(Opt$par)]
    
    
    
    ### Update tau_S, tau_R, tau_U, and tau_V
    S1 = S1 - mean(S1)
    R1 = R1 - mean(R1)
    S2 = S2 - mean(S2)
    R2 = R2 - mean(R2)
    #
    SR1 = cbind(S1,R1)
    SR2 = cbind(S2,R2)
    # 
    UV = cbind(U,V)
    
    
    sigma_UVI = (diag(p) %x% sigma_UV+ matrix(rowSums(apply(UV,1,fun_helper)),
                                              ncol = 2*p, nrow = 2*p))/(upsilon_UV + n - 1 - 2*p)
    
    sigma_UV = sigma_UVI[1:2,1:2]
    sigma_SR1 = (Phi_SR1 + matrix(rowSums(apply(SR1,1,fun_helper)),
                                  ncol = 2, nrow = 2))/(upsilon_SR1 + n - 3)
    sigma_SR2 = (Phi_SR2 + matrix(rowSums(apply(SR2,1,fun_helper)),
                                  ncol = 2, nrow = 2))/(upsilon_SR2 + n - 3)
    
    
    
    
    ### Update alpha
    alpha <- exp(optimize(function(l){
      a <- exp(l)
      obj <- (a_a - 1)*log(a) - b_a*a + lgamma(K*a) - K*lgamma(a) +
        sum((a + Pk - 1)*Elog_t_k)
      return(obj)
    }, c(-10,10), maximum = T)$maximum)
    
    
    alphaTld = alpha + colSums(Pmk)
    
    
    
    ### Update lam
    lam = (B_0 + 0.5*sum(apply(Lambda, 3,
                               function(x) norm(diag(x),type = "2"))))/(A_0 + K*p/2 + 1)
    
    
    setTxtProgressBar(pb,it)
  }
  
  
  ######################################################################################
  ###
  ###       END of INTERATION
  ###
  #######################################################################################
  
  
  Z_est = apply(Pmk_all,1,find_clust)
  # (ARI_weighted = adjustedRandIndex(Z_est, E(A)$cluster))
  # (ARI_init = adjustedRandIndex(Z_init, E(A)$cluster))
  # Z_temp <- ifelse(E(A)$cluster == 5, 1, 2)
  # Z_est2 <- ifelse(Z_est == 1, 1, 2)
  # adjustedRandIndex(Z_est2, Z_temp)
  
  
  ret = list(estimates = list(S1=S1,
                              R1=R1,
                              S2=S2,
                              R2=R2,
                              U = U,
                              V = V,
                              Y = Y,
                              Eta = Eta,
                              Pmk = Pmk_all,
                              alpha = alpha,
                              Elog_t_k = Elog_t_k,
                              E_t_0 = E_t_0,
                              phi = phi,
                              sigma_SR1 = sigma_SR1,
                              sigma_SR2 = sigma_SR2,
                              sigma_UV = sigma_UV,
                              Lambda = Lambda,
                              Z_est = Z_est,
                              Z_init = Z_init),
             userInputs = list(K = K,
                               p = p,
                               Phi_SR1 = Phi_SR1,
                               Phi_SR2 = Phi_SR2,
                               Phi_UV = Phi_UV,
                               upsilon_SR1 = upsilon_SR1,
                               upsilon_SR2 = upsilon_SR2,
                               upsilon_UV = upsilon_UV,
                               a_0 = a_0,
                               A_0 =  A_0,
                               B_0 = B_0,
                               nu_0 = nu_0,
                               eta_0 = eta_0,
                               lambda_a = lambda_a))
  class(ret) = "eClust_GEM"
  return(ret)
  
  
  
  
  
}


