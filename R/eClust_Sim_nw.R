
#' @import igraph
#' @import actuar
#' @import LaplacesDemon
#' @import SimDesign
#' @import dplyr
#' @import truncnorm
#' @export




if (FALSE){
  library(igraph)
  library(actuar)
  library(SimDesign)
  library(LaplacesDemon)
  library(tibble)
  library(truncnorm)
  
  Set=1
  nSims = 1
  n=400
  M=2e4
  showPB=TRUE
  distribution = "Normal"
  Lambda_sc = 5
  set.seed(1234)
  p=2
  K=4
  UVShape=100
  UVRate=40
  mag=0.5
  Ymag = 1
  if_weight = TRUE
  with_noise = TRUE
  noise_mean = 0.5
  noise_prop = 0.1
  phi_min = 0.05
  phi_max = 0.5
}





eClust_Sim_nw = function(Set=1, nSims = 1,n=400,M=2e4,showPB=TRUE,p = 2,
                        K = 2*p,distribution = "Normal",Lambda_sc = 5,
                        UVShape=100, UVRate=40, mag=0.5,Ymag = 1, if_weight = TRUE,
                        with_noise = TRUE, noise_prop = 0.1, noise_mean = 0.5,
                        phi_min = 0.05, phi_max = 0.5){
  
  ZSims = ASims = list()
  nVec = 1:n
  
  if(Set == 1){
    conc = 50
    # UVShape=100;UVRate=40
    Y = c(0,1/2,1,3/2)*pi; K = length(Y)
    UDir = sample(Y,size=n,replace=T);UDir=cbind(cos(UDir),sin(UDir))
    tauSR1=0.5
    tauSR2=0.5
    Y = cbind(cos(Y),sin(Y))*Ymag
  }
  if(Set == 2){
    conc = 50
    # UVShape=20;UVRate=5
    Y = c(0,3/8,9/8,3/2)*pi; K = length(Y)
    UDir = sample(c(3/2,0,3/4)*pi,size=n,replace=T);UDir=cbind(cos(UDir),sin(UDir))
    tauSR1=0.5
    tauSR2=0.5
    Y=cbind(cos(Y),sin(Y))*Ymag
  }
  if(Set == 3){
    conc = 50
    # UVShape=20;UVRate=4
    Y = c(0,2/3,4/3)*pi; K = length(Y)
    UDir = sample(c(1/3,1,5/3)*pi,size=n,replace=T);UDir=cbind(cos(UDir),sin(UDir))
    tauSR1=0.5
    tauSR2=0.5
    Y=cbind(cos(Y),sin(Y))*Ymag
  }
  if(Set == 4){
    conc = 0
    # UVShape=10;UVRate=4
    Y = c(0,1/2,1,3/2)*pi; K = length(Y)
    UDir = runif(n,0,2*pi);UDir=cbind(cos(UDir),sin(UDir))
    tauSR1=0.5
    tauSR2=0.5
    Y = cbind(cos(Y),sin(Y))*Ymag
  }
  if(Set == 5){
    conc = 50
    # UVShape=20;UVRate=4
    Y = rbind(c(1,0,0),
              c(0,1,0),
              c(-1,0,0),
              c(0,-1,0),
              c(0,0,1),
              c(0,0,-1))
    K = nrow(Y)
    UDir = Y[sample(K,size=n,replace=T),]
    tauSR1=0.5
    tauSR2=0.5
    Y=Y*Ymag
  }
  if(Set == 6){
    conc = 0
    # UVShape=20;UVRate=4
    Y = rbind(c(1,0,0),
              c(0,1,0),
              c(-1,0,0),
              c(0,-1,0),
              c(0,0,1),
              c(0,0,-1))
    K = nrow(Y)
    UDir = movMF::rmovMF(n,rep(0,ncol(Y)))
    tauSR1=0.5
    tauSR2=0.5
    Y=Y*Ymag
  }
  if(Set == 7){
    conc = 50
    # UVShape=15;UVRate=4
    Y = c(0,1/2,1,3/2)*pi; K = length(Y)
    UDir = sample(Y,size=n,replace=T);UDir=cbind(cos(UDir),sin(UDir))
    tauSR1=0.5
    tauSR2=0.5
    Y = cbind(cos(Y),sin(Y))*Ymag
  }
  
  # if(Set < 5) Y = cbind(cos(Y),sin(Y))*0.8
  # plot(W,asp=1);plotrix::draw.circle(0,0,1)
  
  Alpha = rep(1/K,K)
  if(showPB) pb = txtProgressBar(0,nSims,style=3)
  for(it in 1:nSims){
    U = V = matrix(0,n,ncol(Y))
    UVMagnitude = rgamma(n,shape=UVShape,rate=UVRate)
    for(i in 1:n){
      UVTheta = movMF::rmovMF(2,theta=conc*UDir[i,])
      U[i,] = UVMagnitude[i] * UVTheta[1,]
      V[i,] = UVMagnitude[i] * UVTheta[2,]
    }
    
    SR1 = rmvnorm(n,sigma=tauSR1*outer(1:2,1:2,function(x,y)0.75^abs(x-y)))
    SR2 = rmvnorm(n,sigma=tauSR2*outer(1:2,1:2,function(x,y)0.75^abs(x-y)))
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
    
   beta = 1.5*sample(c(1, 2.5, 4, 5.5),4)
    ## beta <- 3*c(2, 1, -2, -1)
    EL = matrix(NA,M,2)
    
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
    
    if (with_noise){
      for (m in 1:M){
        if(Z[m] == (K+1)){
          EL[m,1] = sample(n, 1, prob = rep(1/n, n))
          EL[m,2] = sample(nVec[-EL[m,1]], 1, prob = rep(1/(n-1), (n-1)) )
          
          W[m] = rexp(1, rate = 1/noise_mean)
          
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
   
    ## pairs(test3)
    
    EL1 = EL[!duplicated(EL),]
    Z = Z[!duplicated(EL)]
    ASims[[it]] = graph_from_edgelist(EL1,directed=TRUE)
    E(ASims[[it]])$cluster = Z
    
    W = W[!duplicated(EL)]
    E(ASims[[it]])$weight = W
    
    if (!if_weight){
      E(ASims[[it]])$weight = rep(1, length(W))
    }
    
    if(showPB) setTxtProgressBar(pb,it)
  }
  
  return(list(A = ASims, Z = Z, S1=SR1[,1], R1= SR1[,2], U=U, V=V, Y=Y, S2 = SR2[,1], R2 = SR2[,2],
              lam = lam, Lambda = Lambda, beta = beta,phi = phi, Eta = Eta))
}

# if (FALSE){
#   library(dplyr)
#   library(MetBrewer)
#   library(RColorBrewer)
#   
#   set.seed(1234)
#   
#   network_data_1 <- eClust_Sim_nw(Set=1, distribution = "Normal", Lambda_sc = 5, 
#                                  UVShape=100, UVRate=40, mag=0.6, Ymag=0.5, noise_mean = 0.1,
#                                  phi_max = 0.2, noise_prop = 0.1)
#   
#   
#   par(mfrow=c(1,2))
#   # network_data_1 <- eClust_Sim(Set=1)
#   
#   table(network_data_1$Z)
#   summary(E(network_data_1$A[[1]])$weight)
#   
#   
#   hist(log(E(network_data_1$A[[1]])$weight), xlab = "Weight")
#   
#   temp <- cbind(log(E(network_data_1$A[[1]])$weight), E(network_data_1$A[[1]])$cluster)
#   colnames(temp) <- c("weight","cluster")
#   temp <- as_tibble(temp)
#   
#   mean_ord <- temp %>% group_by(cluster) %>% 
#     summarise(mean_weight=mean(weight)) %>% 
#     arrange(mean_weight)
#   
#   temp <- temp %>% arrange(match(cluster, mean_ord$cluster))
#   
#   plot(temp$weight, col=met.brewer("Austria")[temp$cluster], pch=temp$cluster, cex=1.5,
#        xlab="", ylab="weight")
#   
#   temp %>% group_by(cluster) %>% summarize(Mean = mean(weight), Vari = sd(weight))
#   
#   saveRDS(network_data_1, "C:/Project/Thesis/temp/Data/noisy_network_sample1.RDS")
#   
#   
#   
#   
#   
#   
#   plot_lnorm = function(m = c(8.7,2.1,-2.1,-6),
#                         v = c(3.6,3.6,3.6,3.6)){
#     
#     
#     
#     xl = qlnorm(0.99, max(m),v[which.max(m)])
#     yl = dlnorm(exp(min(m) - v[which.min(m)]^2), min(m),v[which.min(m)])
#     
#     cols = gray(seq(0.2,0.8,l = length(m)))
#     
#     curve(dlnorm(x,m[1],v[1]),
#           from = 0,
#           to = xl,
#           xlab = '',
#           ylab = '',
#           bty = 'l',
#           lwd = 3,
#           col = cols[1],
#           ylim = c(0,yl))
#     
#     if(length(m) > 1){
#        for(j in 2:length(m)){
#         curve(dlnorm(x,m[j],v[j]),
#               add = T,
#               col = cols[j],
#               lwd = 3)
#         
#       }
#       
#     }
#       
#     curve(dexp(x,1),
#           add = T,
#           lty = 2,
#           col = 2)
#     
#   }
#   
#   
#   
#   
#   
#   plot_norm = function(m = c(8.7,2.1,-2.1,-6),
#                        v = c(3.6,3.6,3.6,3.6)){
#     
#     
#     
#     xl = qnorm(0.99, max(m),v[which.max(m)])
#     yl = dnorm(m[which.min(v)], m[which.min(v)],min(v))
#     
#     cols = gray(seq(0.2,0.8,l = length(m)))
#     
#     curve(dnorm(x,m[1],v[1]),
#           from = -xl,
#           to = xl,
#           xlab = '',
#           ylab = '',
#           bty = 'l',
#           lwd = 3,
#           col = cols[1],
#           ylim = c(0,yl))
#     
#     if(length(m) > 1){
#       for(j in 2:length(m)){
#         curve(dnorm(x,m[j],v[j]),
#               add = T,
#               col = cols[j])
#         
#       }
#       
#     }
#   
#   }
#   
#   
# plot_lnorm(m = c(3.73, 8.22, 1.49, 5.93),
#              v = rep(1.3,4))
#   
# plot_norm(m =c(3.73, 8.22, 1.49, 5.93),
#             v = rep(1.05,4))
#   
#   qexp(0.95,1);qlnorm(0.05,2,0.2)
#   
#   hist(rlnorm(1e4,8,1))
#   
#   
#   
#   
#   
#   
# }

