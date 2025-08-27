// draw_coef.cpp
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec draw_coef_rc(const arma::mat& RHS_decision_c, 
                         const arma::vec& y_decision_c, 
                         const arma::mat& prior_mat_gamma, 
                         double sigma_true) {
  // Compute (X^T X)/sigma_true + prior_mat_gamma
  arma::mat RHS_prod_decision = (RHS_decision_c.t() * RHS_decision_c) / sigma_true + prior_mat_gamma;
  
  // Compute gamma_mean_c = solve(RHS_prod_decision, X^T y) / sigma_true
  arma::vec gamma_mean_c = arma::solve(RHS_prod_decision, (RHS_decision_c.t() * y_decision_c)) / sigma_true;
  
  // Obtain the Cholesky decomposition: returns an upper triangular matrix U such that U'U = RHS_prod_decision.
  arma::mat chol_factor = arma::chol(RHS_prod_decision);
  
  // Number of coefficients
  int k = RHS_prod_decision.n_cols;
  
  // Draw standard normal random vector of length k
  arma::vec z = arma::randn<arma::vec>(k);
  
  // Solve the triangular system: U * x = z to mimic backsolve.
  // arma::solve with arma::trimatu ensures it treats U as upper triangular.
  arma::vec noise = arma::solve(arma::trimatu(chol_factor), z);
  
  // Combine with the mean draw to get the posterior coefficient draw
  arma::vec gamma_draw_c = gamma_mean_c + noise;
  
  return gamma_draw_c;
}
