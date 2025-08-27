// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' Posterior mean & covariance for a single-outcome regression
//'
//' @param y           Numeric vector (m)     – latent responses
//' @param Z           Numeric matrix (m × k) – full design matrix (already cbind(X , F))
//' @param sigma2      Error variance σ²
//' @param lambda_diag Numeric vector (k)     – diagonal of prior precision Λ₀
//' @param lambda_mu   (optional) Numeric vector (k) – Λ₀·μ₀.  
//'                    If missing or NULL, a 0-vector is assumed.
//'
//' @return A list with components
//'   \item{mu}{posterior mean  (k)  as arma::vec}
//'   \item{V}{posterior covariance (k × k)  as arma::mat}
//'
//' @examples
//' set.seed(1)
//' m <- 100; p <- 3
//' Z <- matrix(rnorm(m * p), m, p)
//' y <- Z %*% rnorm(p) + rnorm(m)
//' lambda_diag <- rep(1/10, p)
//' res <- posterior_mean(y, Z, 1, lambda_diag)  # lambda_mu omitted
//'
//' @export
// [[Rcpp::export]]
Rcpp::List posterior_mean(const arma::vec&            y,
                          const arma::mat&            Z,
                          const double                sigma2,
                          const arma::vec&            lambda_diag,
                          Rcpp::Nullable<arma::vec>   lambda_mu = R_NilValue)
{
    const std::size_t m = y.n_elem;
    const std::size_t k = Z.n_cols;

    if (Z.n_rows != m)
        Rcpp::stop("y and Z must have the same number of rows.");
    if (lambda_diag.n_elem != k)
        Rcpp::stop("lambda_diag must have length equal to ncol(Z).");

    arma::vec lambdaMu;
    if (lambda_mu.isNull()) {
        lambdaMu.zeros(k);                      // default: zero vector
    } else {
        lambdaMu = Rcpp::as<arma::vec>(lambda_mu);
        if (lambdaMu.n_elem != k)
            Rcpp::stop("lambda_mu must have length equal to ncol(Z).");
    }

    const double eps = 1.0 / sigma2;            // precision

    // -------- sufficient statistics -----------------------------
    arma::vec mu_tilde = eps * (Z.t() * y) + lambdaMu;
    arma::mat Vinv     = eps * (Z.t() * Z);
    Vinv.diag()       += lambda_diag;

    // -------- posterior covariance V = Vinv⁻¹ -------------------
    arma::mat V = arma::inv_sympd(Vinv);

    // -------- posterior mean μ ----------------------------------
    arma::vec mu = V * mu_tilde;

    return Rcpp::List::create(Rcpp::Named("mu") = mu,
                              Rcpp::Named("V")  = V);
}
