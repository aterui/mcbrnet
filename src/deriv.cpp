#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

// [[Rcpp::export("deriv")]]
SEXP deriv(double t, arma::vec n, Rcpp::List parms) {

  // R compatible
  // deriv <- function(t, n, parms) {
  //   with(parms, {
  //     dn <- n * (r + A %*% n) + C %*% n
  //     list(dn)
  //   })
  // }

  // parameters
  // r: [vector] ns x np vector for intrinsic growth
  // A: [matrix] (ns x np) x (ns x np) interaction matrix
  // C: [matrix] (ns x np) x (ns x np) connectivity matrix
  // psi: [scalar] disturbance incidence
  // input: [function] approximation function
  arma::vec r = parms["r"];
  arma::mat A = parms["A"];
  arma::mat C = parms["C"];

  // output
  arma::vec dn(n.size());

  // model
  dn = n % (r + A * n) + C * n;

  return(Rcpp::List::create(dn));
}
