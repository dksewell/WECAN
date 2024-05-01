#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::colvec dQSR2(const arma::colvec & SS1,
                   const arma::colvec & SS2,
                   const arma::mat & U,
                   const arma::colvec & RR1,
                   const arma::colvec & RR2,
                   const arma::mat & V,
                   const arma::mat & Y,
                   const arma::cube & Lambda,
                   const arma::colvec & beta,
                   const arma::colvec & phi,
                   const double & lam,
                   const arma::colvec & WW,
                   const arma::mat & Pmk,
                   const arma::mat & Pmki1,
                   const arma::mat & Pmki2,
                   const arma::rowvec & Pk,
                   const IntegerVector & Mi1Index,
                   const IntegerVector & Mi2Index,
                   const IntegerMatrix & Mi1,
                   const IntegerMatrix & Mi2,
                   const IntegerMatrix & EE,
                   const arma::mat & Eta,
                   const arma::mat & inv_sigma_SR1,
                   const arma::mat & inv_sigma_SR2,
                   const arma::mat & inv_sigma_UV,
                   const double & nu_0,
                   const double & eta_0){
  int M = EE.nrow();
  int K = Y.n_rows;
  int p = U.n_cols;
  int n = U.n_rows;
  
  arma::mat A_eta = Eta % Eta/2;
  arma::mat A_eta_dev = Eta;
  arma::mat Pr_eta_dev = arma::zeros(M,K);
  arma::colvec a_phi = phi % phi;
  arma::colvec a_phi_dev = 2*phi;
  arma::colvec pi_phi_dev = - (nu_0 + 1) * phi/(nu_0*eta_0*eta_0 + phi % phi);
  arma::mat Pr_w_eta = arma::ones(M,K);
  arma::colvec pi_phi = pow((1 + phi % phi/(nu_0*eta_0*eta_0)), (0.5*(nu_0 + 1)));
  
  arma::mat W_adj_exp = WW * arma::ones(1,K);
  arma::mat log_h_phi = - W_adj_exp % W_adj_exp * diagmat(1/(phi % phi)) * 0.5 - 
    arma::ones(M, K)* diagmat(log(phi));
  arma::mat h_phi_dev = W_adj_exp % W_adj_exp * diagmat(1/(phi % phi % phi)) - 
    arma::ones(M, K)* diagmat(1/phi);
  
  arma::mat eUiWk = exp(SS1*arma::ones(1,K) + U*Y.t());
  arma::mat eViWk = exp(RR1*arma::ones(1,K) + V*Y.t());
  arma::rowvec fuk = sum(eUiWk,0);
  arma::rowvec fvk = sum(eViWk,0);
  arma::colvec hk = arma::zeros(K);
  
  for(int m=0; m<M; m++){
    for(int k=0; k<K; k++){
      hk(k) = hk(k) + 
        Pmk(m,k)/(fuk(k) - eViWk(EE(m,0) - 1,k));
    }
  }
  
  arma::mat Suk = arma::zeros(K,p);
  arma::mat Svk = arma::zeros(K,p);
  for(int i=0;i<n;i++){
    Suk = Suk + eUiWk.row(i).t()*U.row(i);
    Svk = Svk + eViWk.row(i).t()*V.row(i);
  }
  
  
  arma::mat gradS2 = arma::zeros(n,1);
  
  arma::mat gradR2 = arma::zeros(n,1);
 
  arma::cube gradLambda = arma::zeros(p,p,K);
  
  
  
  for(int i=0;i<n;i++){
    for(int k=0; k<K; k++){
      
      if(Mi1Index(i)>0){
        if(Mi1Index(i) == 1){
          
          gradS2.row(i) = gradS2.row(i) + Pmk(Mi1(i,0) - 1, k)*((WW(Mi1(i,0) - 1)- 
            A_eta_dev(Mi1(i,0) - 1, k))/a_phi(k) - Pr_eta_dev(Mi1(i,0) - 1, k) );
        
          
        }else{
          for(int j=0;j<Mi1Index(i);j++){
            
            gradS2.row(i) = gradS2.row(i) + Pmk(Mi1(i,j) - 1, k)*((WW(Mi1(i,j) - 1)- 
              A_eta_dev(Mi1(i,j) - 1, k))/a_phi(k) - Pr_eta_dev(Mi1(i,j) - 1, k) );
            
           
          }
        }
      }
      
     
      if(Mi2Index(i)>0){
        if(Mi2Index(i) == 1){
          
          gradR2.row(i) = gradR2.row(i) + Pmk(Mi2(i,0) - 1, k)*((WW(Mi2(i,0) - 1) -
            A_eta_dev(Mi2(i,0) - 1, k))/a_phi(k)  - Pr_eta_dev(Mi2(i,0) - 1, k) );
         
          
          
        }else{
          for(int j=0;j<Mi2Index(i);j++){
            gradR2.row(i) = gradR2.row(i) + Pmk(Mi2(i,j) - 1, k)*((WW(Mi2(i,j) - 1) -
              A_eta_dev(Mi2(i,j) - 1, k))/a_phi(k)  - Pr_eta_dev(Mi2(i,j) - 1, k) );
            
          }
        }
      }
      
    
    }
    
    
    gradS2.row(i) = gradS2.row(i) - inv_sigma_SR2(0,0)*SS2(i) - inv_sigma_SR2(0,1)*RR2(i);
   
   
    gradR2.row(i) = gradR2.row(i)- inv_sigma_SR2(1,1)*RR2(i) - inv_sigma_SR2(0,1)*SS2(i);
    
    
  }
  
  for(int k=0; k<K; k++){
    gradLambda.slice(k) = - Lambda.slice(k)/lam;
    
    for(int m=0; m<M; m++){
      
      gradLambda.slice(k) += Pmk(m,k)*((WW(m) -
        A_eta_dev(m, k))/a_phi(k)*diagmat(U.row(EE(m,0) - 1))*diagmat(V.row(EE(m,1) - 1)) -
        Pr_eta_dev(m, k)*diagmat(U.row(EE(m,0) - 1))*diagmat(V.row(EE(m,1) - 1)));
      
  
    }
  }
  
  arma::colvec Lambdavector;
  for(int k=0; k<K; k++){
    Lambdavector = join_cols(Lambdavector, gradLambda.slice(k).as_col());
  }
  
  
  arma::colvec results = join_cols(gradS2.as_col(), gradR2.as_col());
  results = join_cols(results, Lambdavector);

  return results;
}
