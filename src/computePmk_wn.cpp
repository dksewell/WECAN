#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat computePmk_wn(const arma::colvec & SS1,
                       const arma::colvec & RR1,
                       const arma::mat & U,
                       const arma::mat & V,
                       const arma::mat & Y,
                       const IntegerMatrix & EE,
                       const arma::colvec & WW,
                       const arma::mat & Eta,
                       const arma::mat & A_eta,
                       const arma::colvec & a_phi,
                       const arma::mat & Pr_w_eta,
                       const arma::mat & h_phi,
                       const double & lambda_a,
                       const arma::rowvec & Elog_t_k,
                       const double & E_t_0){
  int M = EE.nrow();
  int K = Y.n_rows;
  int n = SS1.n_elem;
  
  arma::mat eUiWk = exp(SS1*arma::ones(1,K) + U*Y.t());
  arma::mat eViWk = exp(RR1*arma::ones(1,K) + V*Y.t());
  
  arma::rowvec fuk = sum(eUiWk,0);
  arma::rowvec fvk = sum(eViWk,0);
  
  arma::mat Pmk = arma::zeros(M,(K+1));
  
  
  for(int m = 0;m<M;m++){
    Pmk(m, 0) = E_t_0/lambda_a/(n*(n-1))*exp(-lambda_a*WW(m));
    
    Pmk.submat(m, 1, m, K) = (1-E_t_0)*h_phi.row(m)%(exp((Eta.row(m)*WW(m) - A_eta.row(m))/ a_phi.t()) %
      eUiWk.row(EE(m,0) - 1)%eViWk.row(EE(m,1) - 1)) /
      fuk/(fvk - eViWk.row(EE(m,0) - 1))/Pr_w_eta.row(m)%exp(Elog_t_k);
    
    Pmk.row(m) = Pmk.row(m)/sum(Pmk.row(m));
  }
  
  return Pmk;
}

