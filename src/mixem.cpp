#include "mixem.h"
#include "misc.h"

using namespace Rcpp;
using namespace arma;

// FUNCTION DEFINITIONS
// --------------------
// Perform several EM updates of the mixture weights.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List mixem_rcpp (const arma::mat& L, const arma::vec& w, const arma::vec& x0,
		 int numiter, double zerothresholdsolution, bool verbose) {
  vec obj(numiter);
  vec nnz(numiter);
  vec dmax(numiter);
  int m = L.n_cols;
  mat P = L;
  vec x = x0;
  vec xold(m);
  vec d(m);
  for (unsigned int i = 0; i < numiter; i++) {
    xold = x;
    
    // Run a single EM update.
    mixem_update(L,w,x,P);
    
    // Record the algorithm's progress.
    nnz(i)  = sum(x > zerothresholdsolution);
    d       = abs(x - xold);
    dmax(i) = d.max();
    if (verbose) {
      // TO DO.
    }
  }
  
  return List::create(Named("x")         = x,
		      Named("objective") = obj,
		      Named("nnz")       = nnz,
		      Named("max.diff")  = dmax);
}


// Perform a single EM update of the mixture weights.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec mixem_update_rcpp (const arma::mat& L, const arma::vec& w,
			     const arma::vec& x0) {
  mat P = L;
  vec x = x0;
  mixem_update(L,w,x,P);
  return x;
}

// Perform a single EM update of the mixture weights.
void mixem_update (const mat& L, const vec& w, vec& x, mat& P) {
  double e = 1e-15;
  
  // Compute the n x m matrix of posterior mixture assignment
  // probabilities (L is an n x m matrix). This is the "E step".
  //
  // The equivalent R code when e = 0 is
  //
  //   P <- t(t(L) * x)
  //   P <- P / rowSums(P)
  // 
  P = L;
  scalecols(P,x);
  normalizerowsbymax(P);
  P += e;
  normalizerows(P);
    
  // Update the mixture weights. This is the "M step".
  x = trans(P) * w;
}
