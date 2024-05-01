#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
double evalMargPost_w(const arma::colvec & SS1,
                    const arma::colvec & RR1,
                    const arma::mat & SR1,
                    const arma::mat & SR2,
                    const arma::mat & UV,
                    const arma::mat & U,
                    const arma::mat & V,
                    const arma::mat & Y,
                    const arma::mat & Eta,
                    const arma::mat & A_eta,
                    const arma::colvec & a_phi,
                    const arma::colvec & pi_phi,
                    const arma::mat & Pr_w_eta,
                    const arma::mat & log_h_phi,
                    const arma::mat & sigma_SR1,
                    const arma::mat & sigma_SR2,
                    const arma::mat & sigma_UV,
                    const arma::rowvec & alph,
                    const IntegerMatrix & EE,
                    const double & upsilon_SR1,
                    const double & upsilon_SR2,
                    const double & upsilon_UV,
                    const arma::mat & Phi_SR1,
                    const arma::mat & Phi_SR2,
                    const arma::mat & Phi_UV,
                    const arma::colvec & WW,
                    const double & a_0,
                    const arma::cube & Lambda,
                    const double & lam,
                    const double & A_0,
                    const double & B_0,
                    const arma::mat & Pmk){
  int M = EE.nrow();
  int K = Y.n_rows;
  int n = U.n_rows;
  int p = U.n_cols;
  
  arma::mat UiWk = SS1*arma::ones(1,K) + U*Y.t();
  arma::mat ViWk = RR1*arma::ones(1,K) + V*Y.t();
  arma::mat eUiWk = exp(UiWk);
  arma::mat eViWk = exp(ViWk);
  arma::rowvec fuk = sum(eUiWk,0);
  arma::rowvec fvk = sum(eViWk,0);
  
  double ret = 0;
  
  double sum_lambda = 0;
  double sum_p_alpha = 0;
  double sum_log_phi = 0;
  
  for(int k = 0; k<K; k++){
  for(int m = 0;m<M;m++){
    ret = ret + Pmk(m,k)*(log_h_phi(m,k) + (Eta(m,k)*WW(m) - A_eta(m,k))/ a_phi(k) - 
      log(Pr_w_eta(m,k)) + 
      UiWk((EE(m,0) - 1),k) + 
      ViWk((EE(m,1) - 1),k) - 
      log(fuk(k)) - log(fvk(k) - eViWk((EE(m,0) - 1),k)) );
  }
  
    sum_p_alpha = sum_p_alpha + (A_0 + sum(Pmk.col(k)) - 1)*log(alph(k)) ;
    sum_lambda = sum_lambda + accu(Lambda.slice(k)%Lambda.slice(k));
    sum_log_phi = sum_log_phi + log(pi_phi(k));
  }

  
  
  ret = ret -
    0.5*trace(SR1*inv(sigma_SR1)*SR1.t()) - 0.5*trace(SR2*inv(sigma_SR2)*SR2.t()) -
    0.5*trace(UV*inv(kron(diagmat(ones(1,p)),sigma_UV))*UV.t()) -
    0.5*(n + p + upsilon_SR1 + 1)*log(det(sigma_SR1)) -
    0.5*(n + p + upsilon_SR2 + 1)*log(det(sigma_SR2)) -
    0.5*(n*p + p + upsilon_UV + 1)*log(det(sigma_UV)) -
    0.5*accu(Y%Y) - 0.5*sum_lambda/lam - B_0/lam -
    (A_0 + 1 + K*p*0.5)*log(lam) +
    sum_p_alpha - 0.5*trace(Phi_SR1*inv(sigma_SR1)) -
    0.5*trace(Phi_SR2*inv(sigma_SR2)) - 0.5*trace(Phi_UV*inv(sigma_UV)) + sum_log_phi;

  return ret;
}