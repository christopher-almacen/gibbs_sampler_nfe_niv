// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;

//' Compute the Cross-Product t(A) * A using Armadillo.
//'
//' @param A A numeric matrix.
//' @return The cross-product of A.
//' @export
// [[Rcpp::export]]
arma::mat crossprod_armadillo_single(const arma::mat &A) {
  // Computes t(A) * A
  return A.t() * A;
}

//' Compute the Cross-Product t(A) * B using Armadillo.
//'
//' @param A A numeric matrix.
//' @param B A numeric matrix.
//' @return The cross-product of A and B.
//' @export
// [[Rcpp::export]]
arma::mat crossprod_armadillo_double(const arma::mat &A, const arma::mat &B) {
  // Computes t(A) * B
  return A.t() * B;
}
