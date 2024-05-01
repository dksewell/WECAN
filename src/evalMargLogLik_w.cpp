#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
double evalMargLogLik_w(const arma::colvec & SS1,
                        const arma::colvec & RR1,
                        const arma::mat & U,
                        const arma::mat & V,
                        const arma::mat & Y,
                        const arma::mat & Eta,
                        const arma::colvec & WW,
                        const arma::rowvec & alph,
                        const arma::colvec & phi,
                        const IntegerMatrix & EE,
                        const arma::mat & Pmk){
  int M = EE.nrow();
  int K = Y.n_rows;
  
  arma::mat UiWk = SS1*arma::ones(1,K) + U*Y.t();
  arma::mat ViWk = RR1*arma::ones(1,K) + V*Y.t();
  arma::mat eUiWk = exp(UiWk);
  arma::mat eViWk = exp(ViWk);
  arma::rowvec fuk = sum(eUiWk,0);
  arma::rowvec fvk = sum(eViWk,0);
  
  arma::mat W_adj_exp = WW * arma::ones(1,K);
  arma::mat log_h_phi = - W_adj_exp % W_adj_exp * diagmat(1/(phi % phi)) * 0.5 - 
    arma::ones(M, K)* diagmat(log(phi));
  arma::mat A_eta = Eta % Eta/2;
  arma::colvec a_phi = phi % phi;
  arma::mat Pr_w_eta = arma::ones(M,K);
  
  double ret = 0;
  
  for(int k = 0; k<K; k++){
    for(int m = 0;m<M;m++){
      ret = ret + Pmk(m,k)*(log_h_phi(m,k) + (Eta(m,k)*WW(m) - A_eta(m,k))/ a_phi(k) - 
        log(Pr_w_eta(m,k)) + 
        UiWk((EE(m,0) - 1),k) + 
        ViWk((EE(m,1) - 1),k) - 
        log(fuk(k)) - log(fvk(k) - eViWk((EE(m,0) - 1),k)) );
    }
  }
  return ret;
}
