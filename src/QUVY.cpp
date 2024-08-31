#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double QUVY(const arma::colvec & SS1,
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
            const arma::rowvec & Pk,
            const IntegerMatrix & EE,
            const arma::mat & Eta,
            const arma::mat & sigma_SR1,
            const arma::mat & sigma_SR2,
            const arma::mat & sigma_UV,
            const arma::mat & SR1,
            const arma::mat & SR2,
            const arma::mat & UV,
            const double & nu_0,
            const double & eta_0){
  
  int M = EE.nrow();
  int K = Y.n_rows;
  int p = U.n_cols;

  
  arma::mat A_eta = Eta % Eta/2;
  // arma::mat A_eta_dev = Eta;
  arma::mat Pr_eta_dev = arma::zeros(M,K);
  arma::colvec a_phi = phi % phi;
  // arma::colvec a_phi_dev = 2*phi;
  // arma::colvec pi_phi_dev = (nu_0 + 1) * phi/(nu_0*eta_0*eta_0 + phi % phi);
  arma::mat Pr_w_eta = arma::ones(M,K);
  arma::colvec pi_phi = pow(((1 + phi % phi/(nu_0*eta_0*eta_0))), - (0.5*(nu_0 + 1)));
  
  arma::mat W_adj_exp = WW * arma::ones(1,K);
  arma::mat log_h_phi =  - W_adj_exp % W_adj_exp * diagmat(1/(phi % phi)) * 0.5 - 
    arma::ones(M, K)* diagmat(log(phi));
  // arma::mat h_phi_dev = W_adj_exp% W_adj_exp * diagmat(1/(phi % phi % phi)) - arma::ones(M, K)* diagmat(1/phi);
  
  arma::mat eUiWk = exp(SS1*arma::ones(1,K) + U*Y.t());
  arma::mat eViWk = exp(RR1*arma::ones(1,K) + V*Y.t());
  arma::rowvec fuk = sum(eUiWk,0);
  arma::rowvec fvk = sum(eViWk,0);
  
  double ret = 0;
  
  for(int k=0; k<K; k++){
    for(int m = 0;m<M;m++){
      ret += Pmk(m,k)*(log_h_phi(m,k) + (Eta(m,k)*WW(m) - A_eta(m,k))/ a_phi(k) -
        log(Pr_w_eta(m,k)) + SS1(EE(m,0) - 1) + RR1(EE(m,1) - 1) +
      dot((U.row(EE(m,0) - 1) + V.row(EE(m,1) - 1)), Y.row(k))-
      log(fuk(k)) - log( fvk(k) - eViWk(EE(m,0) - 1,k) ) ) ;
    }
    ret += - 0.5*accu(Lambda.slice(k)%Lambda.slice(k))/lam;
    ret += log(pi_phi(k));
  }
  
  ret += - 0.5*trace(SR1*inv(sigma_SR1)*SR1.t()) - 0.5*trace(SR2*inv(sigma_SR2)*SR2.t())-
    0.5*trace(UV*inv(kron(diagmat(ones(1,p)),sigma_UV))*UV.t()) - 0.5*accu(Y%Y);
  
  return ret ;
}
