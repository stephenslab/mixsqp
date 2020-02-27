#ifndef INCLUDE_MISC
#define INCLUDE_MISC

#include <RcppArmadillo.h>

// FUNCTION DECLARATIONS
// ---------------------
void scalecols          (arma::mat& A, const arma::vec& b);
void normalizerows      (arma::mat& A);
void normalizerowsbymax (arma::mat& A);

#endif
