#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double computeELBO_wn(const arma::colvec & SS1,
                      const arma::colvec & RR1,
                      const arma::colvec & beta,
                      const arma::mat & SR1,
                      const arma::mat & SR2,
                      const arma::mat & UV,
                      const arma::mat & U,
                      const arma::mat & V,
                      const arma::mat & Y,
                      const arma::mat & Eta, 
                      const arma::mat & log_h_phi,
                      const arma::mat & A_eta,
                      const arma::colvec & a_phi,
                      const arma::colvec & pi_phi,
                      const arma::mat & Pr_w_eta,
                      const IntegerMatrix & EE,
                      const arma::mat & Phi_SR1,
                      const arma::mat & Phi_SR2,
                      const arma::mat & Phi_UV,
                      const arma::colvec & WW,
                      const double & a_0,
                      const arma::cube & Lambda,
                      const double & lam,
                      const double & A_0,
                      const double & B_0,
                      const arma::mat & Pmk,
                      const arma::rowvec & Pk,
                      const arma::rowvec & Pm0,
                      const arma::rowvec & Elog_t_k,
                      const double & E_t_0,
                      const double & lambda_a,
                      const double & alpha,
                      const double & PmklogPmk,
                      const arma::rowvec & alphaTld,
                      const double x_0){
                        
                      
                      
  
  int M = EE.nrow();
  int K = Y.n_rows;
  int n = SS1.n_elem;
  
  
  arma::mat eUiWk = exp(SS1*arma::ones(1,K) + U*Y.t());
  arma::mat eViWk = exp(RR1*arma::ones(1,K) + V*Y.t());
  arma::rowvec fuk = sum(eUiWk,0);
  arma::rowvec fvk = sum(eViWk,0);
  
  double ret = 0;
  
 
  for(int m = 0;m<M;m++){
    for(int k=0; k<K; k++){
      ret += Pmk(m,k)*(log_h_phi(m,k) + (Eta(m,k)*WW(m) - A_eta(m,k))/ a_phi(k) -
        log(Pr_w_eta(m,k)) + SS1(EE(m,0) - 1) + RR1(EE(m,1) - 1) +
        dot((U.row(EE(m,0) - 1) + V.row(EE(m,1) - 1)), Y.row(k))-
        log(fuk(k)) - log( fvk(k) - eViWk(EE(m,0) - 1,k)) + 
        Elog_t_k(k) + log(1-E_t_0) ) ;
    }
    ret += Pm0(m)*(log(E_t_0) - log(n*(n-1)/lambda_a) - WW(m)*lambda_a);
  }
  
  ret += - x_0*log(E_t_0) - (M-x_0)*log(1-E_t_0) +
    sum((alpha + Pk - alphaTld)%Elog_t_k) -
    PmklogPmk + sum(lgamma(alphaTld));
  
  return ret ;
}

