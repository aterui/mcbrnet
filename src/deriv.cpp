#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

// [[Rcpp::export("deriv")]]
SEXP deriv(double t, arma::vec n, Rcpp::List parms) {

  // R compatible
  // deriv <- function(t, n, parms) {
  //   with(parms, {
  //     psi <- input(t)
  //     dn <- n * (r - psi * e + t(A) %*% n) + C %*% n
  //     list(dn)
  //   })
  // }

  // parameters
  // r: [vector] ns x np vector for intrinsic growth
  // A: [matrix] (ns x np) x (ns x np) interaction matrix
  // C: [matrix] (ns x np) x (ns x np) connectivity matrix
  // e: [vector] ns x np vector for disturbance intensity
  // psi: [scalar] disturbance incidence
  // input: [function] approximation function
  arma::vec r = parms["r"];
  arma::vec e = parms["e"];
  arma::mat A = parms["A"];
  arma::mat C = parms["C"];
  Rcpp::Function input = parms["input"];

  // output
  arma::vec dn(n.size());

  // model
  SEXP v = input(t);
  Rcpp::NumericVector psi(v);
  dn = n % (r - psi[0] * e + A.t() * n) + C * n;

  return(Rcpp::List::create(dn));
}
