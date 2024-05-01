#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat computeEta(const arma::colvec & SS2,
                      const arma::colvec & RR2,
                      const arma::colvec & beta,
                      const arma::mat & U,
                      const arma::mat & V,
                      const IntegerMatrix & EE,
                      const arma::cube & Lambda,
                      const int & K){
  int M = EE.nrow();
  
  arma::mat eta = arma::zeros(M,K);
  arma::mat eta_helper = arma::zeros(1,1);
  
  for(int k = 0;k<K;k++){
    for(int m = 0;m<M;m++){
      eta_helper = SS2(EE(m,0) - 1) + RR2(EE(m,1) - 1) + U.row(EE(m,0) - 1)*Lambda.slice(k)*V.row(EE(m,1) - 1).t();
      eta(m,k) = eta_helper(0,0) + beta(k);
      
    }
  }
  
  return eta;
}