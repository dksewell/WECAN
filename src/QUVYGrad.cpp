#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::colvec dQUVY(const arma::colvec & SS1,
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
  
  
  arma::mat gradS1 = arma::zeros(n,1);
  arma::mat gradS2 = arma::zeros(n,1);
  arma::mat gradU = arma::zeros(n,p);
  
  arma::mat gradR1 = arma::zeros(n,1);
  arma::mat gradR2 = arma::zeros(n,1);
  arma::mat gradV = arma::zeros(n,p);
  arma::colvec gradphi = arma::zeros(K);
  
  arma::mat gradY = -Y;
  arma::colvec gradbeta = arma::zeros(K,1);
  arma::cube gradLambda = arma::zeros(p,p,K);
  
  
  
  for(int i=0;i<n;i++){
    for(int k=0; k<K; k++){
      
      arma::mat VLambda = V*Lambda.slice(k);
      arma::mat ULambda = U*Lambda.slice(k);
 
      
      
      if(Mi1Index(i)>0){
        if(Mi1Index(i) == 1){
          
          gradS2.row(i) = gradS2.row(i) + Pmk(Mi1(i,0) - 1, k)*((WW(Mi1(i,0) - 1)- 
            A_eta_dev(Mi1(i,0) - 1, k))/a_phi(k) - Pr_eta_dev(Mi1(i,0) - 1, k) );
          
          gradU.row(i) = gradU.row(i) + Pmk(Mi1(i,0) - 1, k)*((WW(Mi1(i,0) - 1)- 
            A_eta_dev(Mi1(i,0) - 1, k))/a_phi(k) - 
            Pr_eta_dev(Mi1(i,0) - 1, k) )*VLambda.row(EE(Mi1(i,0) - 1, 1) - 1);
          
        }else{
          for(int j=0;j<Mi1Index(i);j++){
            
            gradS2.row(i) = gradS2.row(i) + Pmk(Mi1(i,j) - 1, k)*((WW(Mi1(i,j) - 1)- 
              A_eta_dev(Mi1(i,j) - 1, k))/a_phi(k) - Pr_eta_dev(Mi1(i,j) - 1, k) );
            
            gradU.row(i) = gradU.row(i) + Pmk(Mi1(i,j) - 1, k)*((WW(Mi1(i,j) - 1)- 
              A_eta_dev(Mi1(i,j) - 1, k))/a_phi(k) - 
              Pr_eta_dev(Mi1(i,j) - 1, k) )*VLambda.row(EE(Mi1(i,j) - 1, 1) - 1);
          }
        }
      }
      
      gradS1.row(i) = gradS1.row(i) + Pmki1(i,k) - Pk(k)*eUiWk(i,k)/fuk(k);
      gradU.row(i) = gradU.row(i) + (Pmki1(i,k) - Pk(k)*eUiWk(i,k)/fuk(k))*Y.row(k);
      
      if(Mi2Index(i)>0){
        if(Mi2Index(i) == 1){

          gradR2.row(i) = gradR2.row(i) + Pmk(Mi2(i,0) - 1, k)*((WW(Mi2(i,0) - 1) -
            A_eta_dev(Mi2(i,0) - 1, k))/a_phi(k)  - Pr_eta_dev(Mi2(i,0) - 1, k) );

          gradV.row(i) = gradV.row(i) + Pmk(Mi2(i,0) - 1, k)*((WW(Mi2(i,0) - 1) -
            A_eta_dev(Mi2(i,0) - 1, k))/a_phi(k) -
            Pr_eta_dev(Mi2(i,0) - 1, k) )*ULambda.row(EE(Mi2(i,0) - 1, 0) - 1);


        }else{
          for(int j=0;j<Mi2Index(i);j++){
            gradR2.row(i) = gradR2.row(i) + Pmk(Mi2(i,j) - 1, k)*((WW(Mi2(i,j) - 1) -
              A_eta_dev(Mi2(i,j) - 1, k))/a_phi(k)  - Pr_eta_dev(Mi2(i,j) - 1, k) );

             gradV.row(i) = gradV.row(i) + Pmk(Mi2(i,j) - 1, k)*((WW(Mi2(i,j) - 1) -
               A_eta_dev(Mi2(i,j) - 1, k))/a_phi(k) -
             Pr_eta_dev(Mi2(i,j) - 1, k) )*ULambda.row(EE(Mi2(i,j) - 1, 0) - 1);
          }
        }
      }
      
      gradR1.row(i) = gradR1.row(i) + Pmki2(i,k) - (eViWk(i,k)*(hk(k) - 
        Pmki1(i,k)/(fvk(k) - eViWk(i,k))));
                                                      
      gradV.row(i) = gradV.row(i) + (Pmki2(i,k) - (eViWk(i,k)*(hk(k) -
        Pmki1(i,k)/(fvk(k) - eViWk(i,k)))))*Y.row(k);
      

    }
    
    gradS1.row(i) = gradS1.row(i) - inv_sigma_SR1(0,0)*SS1(i) - inv_sigma_SR1(0,1)*RR1(i);
    gradS2.row(i) = gradS2.row(i) - inv_sigma_SR2(0,0)*SS2(i) - inv_sigma_SR2(0,1)*RR2(i);
    gradU.row(i) = gradU.row(i) - inv_sigma_UV(0,0)*U.row(i) - inv_sigma_UV(0,1)*V.row(i);
    
    gradR1.row(i) = gradR1.row(i)- inv_sigma_SR1(1,1)*RR1(i) - inv_sigma_SR1(0,1)*SS1(i);
    gradR2.row(i) = gradR2.row(i)- inv_sigma_SR2(1,1)*RR2(i) - inv_sigma_SR2(0,1)*SS2(i);
    gradV.row(i) = gradV.row(i) - inv_sigma_UV(1,1)*V.row(i) - inv_sigma_UV(0,1)*U.row(i);
    

  }
  
  for(int k=0; k<K; k++){
    gradLambda.slice(k) = - Lambda.slice(k)/lam;
    
    for(int m=0; m<M; m++){
  
  gradY.row(k) += Pmk(m,k)*(U.row(EE(m,0) - 1) + V.row(EE(m,1) - 1) -
        Suk.row(k)/fuk(k) -  (Svk.row(k) - eViWk(EE(m,0) - 1,k)*V.row(EE(m,0) - 1))/
          (fvk(k) - eViWk(EE(m,0) - 1,k)));
  
  gradLambda.slice(k) += Pmk(m,k)*((WW(m) -
    A_eta_dev(m, k))/a_phi(k)*diagmat(U.row(EE(m,0) - 1))*diagmat(V.row(EE(m,1) - 1)) -
    Pr_eta_dev(m, k)*diagmat(U.row(EE(m,0) - 1))*diagmat(V.row(EE(m,1) - 1)));
  
    
    gradbeta(k) += Pmk(m,k)*((WW(m) -  A_eta_dev(m, k))/a_phi(k)-Pr_eta_dev(m,k)/Pr_w_eta(m,k));
    
    gradphi(k) += Pmk(m,k)*(h_phi_dev(m,k) - 
      (Eta(m)*WW(m)- A_eta(m))*a_phi_dev(k)/pow(a_phi(k),2) - Pr_eta_dev(m,k));
    }
  }
  
  arma::colvec Lambdavector;
  for(int k=0; k<K; k++){
    Lambdavector = join_cols(Lambdavector, gradLambda.slice(k).as_col());
    gradphi(k) += pi_phi_dev(k);
  }
  
  
  arma::colvec results = join_cols(gradS1.as_col(), gradS2.as_col());
  results = join_cols(results, gradU.as_col());
  results = join_cols(results, gradR1.as_col());
  results = join_cols(results, gradR2.as_col());
  results = join_cols(results, gradV.as_col());
  results = join_cols(results, gradY.as_col());
  results = join_cols(results, Lambdavector);
  
  results = join_cols(results, gradbeta.as_col());
  results = join_cols(results, gradphi.as_col());

  
  return results;
}
