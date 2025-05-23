#' wecan: Weighted Edge Clustering Adjusting for Noise
#' 
#' This function implements the weighted edge clustering algorithm 
#' of Li and Sewell
#' 
#' @param A Network object of class igraph
#' @param K integer. Maximum number of clusters, by default is 20
#' @param p integer. Dimension of latent space, by default is 3
#' @param maxIter integer. Maximum number of iterations of the VB-EM algorithm
#' @param maxIterVB integer. Maximum number of iterations of the 
#' VB step nested within the EM algorithm
#' @param a_0 numeric. Hyper-parameters.
#' @param A_0 numeric. Hyper-parameters.
#' @param B_0 numeric. Hyper-parameters.
#' @param nu_0 numeric. Hyper-parameters.
#' @param eta_0 numeric. Hyper-parameters.
#' @param CGSteps Number of steps in each conjugate gradient update
#' nested within the VB-EM algorithm
#' @param distribution character.  Currently, only "Normal" and "Poisson" are supported
#' @param weight_attr_name character.  Name of the weight attribute of A to be modeled
#' @param eps numeric. Convergence threshold for relative change of ELBO in variational E step
#' @param upsilon_SR1 numeric. Degrees of freedom of the inverse Wishart distribution, which serves as the prior for the covariance matrix of S1 and R1
#' @param upsilon_SR2 numeric. Degrees of freedom of the inverse Wishart distribution, which serves as the prior for the covariance matrix of S2 and R2
#' @param upsilon_UV numeric. Degrees of freedom of the inverse Wishart distribution, which serves as the prior for the covariance matrix of U and V
#' @param Phi_SR1 matrix. The scale matrix of the inverse Wishart distribution, which serves as the prior for the covariance matrix of S1 and R1
#' @param Phi_SR2 matrix. The scale matrix of the inverse Wishart distribution, which serves as the prior for the covariance matrix of S2 and R2
#' @param Phi_UV matrix. The scale matrix of the inverse Wishart distribution, which serves as the prior for the covariance matrix of U and V
#' @param a_a numeric. Hyper-parameters.
#' @param b_a numeric. Hyper-parameters.
#' @param t_0 numeric. The initial proportion of the noisy edges, range between 0 and 1.
#' @param alpha_0 numeric. initial alpha
#' @param beta_0 vector. The initial weights mean for each cluster.
#' 
#' @return Object of class 'wecan' with the following elements:
#' \itemize{
#' \item estimates, a list with the following elements:
#'    \itemize{
#'    \item S1: a vector representing overall propensities of nodes to send a large or small number of edges
#'    \item R1: a vector representing overall propensities of nodes to receive a large or small number of edges
#'    \item S2: a vector representing overall propensities of nodes to send edges with large or small weights
#'    \item R2: a vector representing overall propensities of nodes to receive edges with large or small weights
#'    \item U: a matrix representing nodes' latent sending features
#'    \item V: a matrix representing nodes' latent receiving features
#'    \item Y: a matrix representing edge's cluster features
#'    \item Eta: estimated canonical parameters for the exponential distribution family
#'    \item Pmk: a matrix that describe the probability of each edge belonging to each cluster 
#'    \item phi: estimated dispersion parameter for the exponential distribution family
#'    \item sigma_SR1: covariance matrix for S1 and R1
#'    \item sigma_SR2: covariance matrix for S2 and R2
#'    \item sigma_UV: covariance matrix for U and V
#'    \item Lambda: array representing how each latent cluster features determine the weights
#'    \item Z_est: estimated edge clusters by WECAN model, "1" suggest a noisy cluster
#'    }
#' \item userInputs, a list with the following elements:
#'    \itemize{
#'    \item K Maximum number of clusters.
#'    \item p Dimension of latent space.
#'    \item Phi_SR1 The scale matrix of the inverse Wishart distribution, which serves as the prior for the covariance matrix of S1 and R1
#'    \item Phi_SR2 The scale matrix of the inverse Wishart distribution, which serves as the prior for the covariance matrix of S2 and R2
#'    \item Phi_UV The scale matrix of the inverse Wishart distribution, which serves as the prior for the covariance matrix of U and V
#'    \item upsilon_SR1 Degrees of freedom of the inverse Wishart distribution, which serves as the prior for the covariance matrix of S1 and R1
#'    \item upsilon_SR2 Degrees of freedom of the inverse Wishart distribution, which serves as the prior for the covariance matrix of S2 and R2
#'    \item upsilon_UV Degrees of freedom of the inverse Wishart distribution, which serves as the prior for the covariance matrix of U and V
#'    \item a_0 Hyper-parameters
#'    \item A_0 Hyper-parameters
#'    \item B_0 Hyper-parameters
#'    \item nu_0 Hyper-parameters
#'    \item eta_0 Hyper-parameters
#'    \item lambda_a Hyper-parameters
#'    }
#' \item network, the original igraph object
#' }
#' 
#' @examples
#' \dontrun{
#' #Fill in some code here to show how to use wecan(), BIC(), and ICL().
#' }
#' 
#' @import igraph
#' @import MASS
#' @import LaplacesDemon
#' @import RcppArmadillo
#' @import Rcpp
#' @import aricode
#' @import aLSEC
#' @import mclust
#' @export
#' @exportClass wecan



wecan = function(A, K = 20,
                         p = 3,
                         maxIter=1e3,
                         maxIterVB=100,
                         a_0 = 1,
                         A_0 = 3,
                         B_0 = 0.5,
                         nu_0 = 2,
                         eta_0 = 5,
                         CGSteps=25,
                         distribution = "Normal",
                         weight_attr_name = "weight",
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
  W = igraph::edge_attr(A,weight_attr_name)
  W_adj = W
  A_adj <- as_adjacency_matrix(A, attr = weight_attr_name)
  
  
  ### Index edges
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
  
  
  
  S2 = strength(A, mode = "out")
  S2 = (S2 - mean(S2))/max(abs(S2 - mean(S2)))
  S2_init = S2
  R2 = strength(A, mode = "in")
  R2 = (R2 - mean(R2))/max(abs(R2 - mean(R2)))
  R2_init = R2
  SR2 = cbind(S2,R2)
  
  
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
    
    phi = rep(sd(W_adj), K)
    
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
  
  
  Pm0 = numeric(M)
  Pm0[noise_index] = 1
  Pm0[non_noise_index] = 0
  Pmk = init$estimates$Pmk
  
  
  Pmk_all = matrix(0, nrow = M, ncol = (K+1))
  Pmk_all[,1] =  Pm0
  Pmk_all[non_noise_index,2:(K+1)] = Pmk
  
  Pmk = Pmk_all[,2:(K+1)]
  
  
  Pk = colSums(Pmk)
  alphaTld = alpha + colSums(Pmk)
  x_0 = sum(Pm0)
  
  
  Z_prev <- apply(Pmk_all,1,find_clust)
  
  
  
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
    
    Pmk_all = computePmk_wn(S1,R1,U,V,Y,EE,W_adj,Eta,A_eta,a_phi,Pr_w_eta,h_phi, lambda_a, Elog_t_k, E_t_0)
    Pmk = Pmk_all[,2:(K+1)]
    Pm0 = Pmk_all[,1]
    
    Pk = colSums(Pmk)
    PmklogPmk = sum(Pmk_all[Pmk_all!=0]*log(Pmk_all[Pmk_all!=0]))
    
    Pmki1 = getPmki(Pmk,Mi1Index,Mi1)
    Pmki2 = getPmki(Pmk,Mi2Index,Mi2)
    
    
    alphaTld = alpha + colSums(Pmk)
    x_0 = sum(Pm0)
    
    
    ELBO_Estep = numeric(maxIterVB)
    
    
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
                              Z_est = Z_est),
             userInputs = list(K = K,
                               p = p,
                               weight_attr_name = weight_attr_name,
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
                               lambda_a = lambda_a),
             A = A)
  class(ret) = "wecan"
  
  return(ret)
}


