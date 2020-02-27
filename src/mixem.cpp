#include "mixem.h"
#include "misc.h"

using namespace arma;

// FUNCTION DEFINITIONS
// --------------------
// Perform a single expectation maximization (EM) update.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec mixem_rcpp (const arma::mat& L, const arma::vec& w,
		      const arma::vec& x0) {
  mat P = L;
  vec x = x0;
  mixem_update(L,w,x,P);
  return x;
}

// Perform a single EM update.
void mixem_update (const mat& L, const vec& w, vec& x, mat& P) {
  double e = 1e-15;
  
  // Compute the posterior mixture assignment probabilities. This is
  // the "E step".
  P = L;
  scalecols(P,x);
  normalizerowsbymax(P);
  P += e;
  normalizerows(P);
    
  // Update the mixture weights. This is the "M step".
  x = trans(P) * w;
}
