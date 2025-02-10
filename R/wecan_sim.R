#' Title: Create the simulated networks
#' Description: This function created simulated networks with clustered edges.
#' @param p integer. Number of latent dimensions
#' @param K integer. True number of clusters
#' @param n integer. Number of actors in the network
#' @param M integer. Number of edges in the network
#' @param conc concentration parameter of the mixtures of von Mises-Fisher distributions to generate U, V
#' @param UVShape shape parameter of the Gamma distribution to generate the magnitude of U, V
#' @param UVRate rate parameter of the Gamma distribution to generate the magnitude of U, V
#' @param if_weight logical. indicates whether simulate a weighted network
#' @param distribution character. The distribution for network's weights, currently only "Normal" can be used.
#' @param with_noise logical. indicates whether simulate noisy edges
#' @param noise_mean numeric. mean weights for the noisy edges
#' @param beta vector. mean weights for edges in each cluster
#' @return simulated networks, edge clusters and the nodes-centric/edge-centric features
#' @export 
#'
#' @import movMF
#' @import mvtnorm
#' @import truncnorm

wecan_sim = function(p = 2,
                     K = 2*p,
                     n = 400, Lambda_sc = 5,
                     M = round(n * (n-1) * 0.05),
                     conc = c(1.5, 2.5, 6, 25, 150, 200,
                                       1.75, 5, 10, 50, 150, 200)[6 * (p - 2) + K - 1],
                     UVShape=20,
                     UVRate=4, if_weight = TRUE, mag=1,Ymag = 1,
                     with_noise = TRUE, noise_mean = 0.5,
                     noise_prop = 0.1, distribution = "Normal",
                     phi_min = 0.05, phi_max = 0.5,
                     beta = c(4,0,2,13)){
  
  
  tauSR1=0.5
  tauSR2=0.5
  
  if(p == 2){
    phi = seq(0,2*pi,l = K + 1)[-(K+1)]
    Y = cbind(cos(phi),sin(phi))*Ymag
    rm(phi)
  }else{
    # Create function to generate points roughly equidistant on sphere
    fibonacci = function(N){
      pts = matrix(0.0, N, 3)
      phi = pi * (3 - sqrt(5))
      
      for(i in 0:(N-1)){
        y = 1 - (i / (N-1)) * 2
        radius = sqrt(1 - y^2)
        theta = phi * i
        
        x = cos(theta) * radius
        z = sin(theta) * radius
        
        pts[i+1,] = c(x,y,z)
      }
      
      pts
    }
    
    Y = fibonacci(K)*Ymag
    rm(phi)
    rm(theta)
  }
  
  UDir = Y[sample(K,n,TRUE),]
  U = V = matrix(0,n,p)
  UVMagnitude = rgamma(n,shape=UVShape,rate=UVRate)
  for(i in 1:n){
    UVTheta = movMF::rmovMF(2,theta=conc*UDir[i,])
    U[i,] = UVMagnitude[i] * UVTheta[1,]
    V[i,] = UVMagnitude[i] * UVTheta[2,]
  }
  SR1 = mvtnorm::rmvnorm(n,sigma=tauSR1*outer(1:2,1:2,function(x,y)0.75^abs(x-y)))
  SR2 = mvtnorm::rmvnorm(n,sigma=tauSR2*outer(1:2,1:2,function(x,y)0.75^abs(x-y)))
  
  eUiWk = exp(SR1[,1]%*%matrix(1,1,K) + tcrossprod(U,Y))
  eViWk = exp(SR1[,2]%*%matrix(1,1,K) + tcrossprod(V,Y))
  
  EL = matrix(NA,M,2)
  
  
  lam = 0.5
  # Lambda =replicate(K,diag(rnorm(ncol(U), 0, lam)))
  # beta = 0
  
  if (Lambda_sc==1) { # All Random
    Lambda =replicate(K,diag(rnorm(ncol(U), 0, lam)))
  }
  
  if (Lambda_sc==2){ # All Positive
    Lambda = abs(replicate(K,diag(rnorm(ncol(U), 0, lam))))
  }
  
  if (Lambda_sc == 3){ # All Negative
    Lambda = -abs(replicate(K,diag(rnorm(ncol(U), 0, lam))))
  }
  
  if (Lambda_sc == 4){ #One positive, one negative, others are random
    Lambda = array(NA, dim = c(p,p,K))
    Lambda[,,1] = abs(diag(rnorm(ncol(U), 0, lam)))
    Lambda[,,2] = -abs(diag(rnorm(ncol(U), 0, lam)))
    Lambda[,,3:K] = replicate(K-2,diag(rnorm(ncol(U), 0, lam)))
    
  }
  if (Lambda_sc == 5){ #Given the certain values
    Lambda = array(NA, dim = c(p,p,K))
    Lambda[,,1] = mag*diag(c(1,1))
    Lambda[,,2] = mag*diag(c(-1,1))
    Lambda[,,3] = mag*diag(c(1,-1))
    Lambda[,,4] = mag*diag(c(-1,-1)) 
  }
  
  
  
  
  if (with_noise){
    Pi = c((1-noise_prop)*rep(1/K, K), noise_prop)
    Z = sample( (K+1), size=M, replace=T, prob= Pi)
  }else{
    Alpha = rep(1/K,K)
    Z = sample(K,size=M,replace=T,prob=Alpha)
  }
  
  W = numeric(length(Z))
  Eta = numeric(length(Z))
  phi = runif(K, min = phi_min, max = phi_max)
  
  
  
  nVec = 1:n
  if (with_noise){
    for (m in 1:M){
      if(Z[m] == (K+1)){
        EL[m,1] = sample(n, 1, prob = rep(1/n, n))
        EL[m,2] = sample(nVec[-EL[m,1]], 1, prob = rep(1/(n-1), (n-1)) )
        
        W[m] = rgamma(1, shape = 2, rate = 2/noise_mean)
        
      }else{
        EL[m,1] = sample(n,1,prob=eUiWk[,Z[m]]) 
        EL[m,2] = sample(nVec[-EL[m,1]],1,prob=eViWk[-EL[m,1],Z[m]])
        Eta[m] = SR2[EL[m,1],1] + SR2[EL[m,2],2] + U[EL[m,1],] %*% Lambda[,,Z[m]] %*% V[EL[m,2],] + beta[Z[m]]
        
        if(distribution == "Poisson") {
          phi = rep(0, K)
          W[m] = rztpois(1, exp(Eta[m]))
        }
        if (distribution == "Normal"){
          
          W[m] = rnorm(1, mean = Eta[m], sd = phi[Z[m]])
          W[m] = exp(W[m])
          
        }
        if (distribution == "truncated_normal"){
          W[m] = rtruncnorm(1,a = 0, mean = Eta[m], sd = phi[Z[m]])
        }
      }
    }
  }else{
    for(m in 1:M){
      EL[m,1] = sample(n,1,prob=eUiWk[,Z[m]])
      EL[m,2] = sample(nVec[-EL[m,1]],1,prob=eViWk[-EL[m,1],Z[m]])
      
      Eta[m] = SR2[EL[m,1],1] + SR2[EL[m,2],2] + U[EL[m,1],] %*% Lambda[,,Z[m]] %*% V[EL[m,2],] + beta[Z[m]]
      # test1[m] = SR2[EL[m,1],1] + SR2[EL[m,2],2] 
      # test2[m] = U[EL[m,1],] %*% Lambda[,,Z[m]] %*% V[EL[m,2],]
      # for(k in 1:K)test3[m,k] = U[EL[m,1],] %*% Lambda[,,k] %*% V[EL[m,2],] + beta[Z[m]]
      # 
      if(distribution == "Poisson") {
        phi = rep(0, K)
        W[m] = rztpois(1, exp(Eta[m]))
      }
      if (distribution == "Normal"){
        
        W[m] = rnorm(1, mean = Eta[m], sd = phi[Z[m]])
        W[m] = exp(W[m])
        
      }
      if (distribution == "truncated_normal"){
        W[m] = rtruncnorm(1,a = 0, mean = Eta[m], sd = phi[Z[m]])
      }
      
    }
  }
  
  
  EL1 = EL[!duplicated(EL),]
  Z = Z[!duplicated(EL)]
  ASims = graph_from_edgelist(EL1,directed=TRUE)
  E(ASims)$cluster = Z
  
  W = W[!duplicated(EL)]
  E(ASims)$weight = W
  
  if (!if_weight){
    E(ASims)$weight = rep(1, length(W))
  }
  
  
  return(list(A = ASims, Z = Z, S1=SR1[,1], R1= SR1[,2], U=U, V=V, Y=Y, S2 = SR2[,1], R2 = SR2[,2],
              lam = lam, Lambda = Lambda, beta = beta,phi = phi, Eta = Eta))
}
