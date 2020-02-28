#ifndef INCLUDE_OBJECTIVE
#define INCLUDE_OBJECTIVE

#include <RcppArmadillo.h>

// FUNCTION DECLARATIONS
// ---------------------
double compute_objective (const arma::mat& L, const arma::vec& w,
			  const arma::vec& x, const arma::vec& z,
			  const arma::vec& e);

double compute_objective (const arma::mat& L, const arma::mat& U,
			  const arma::mat& V, const arma::vec& w,
			  const arma::vec& x, const arma::vec& z,
			  const arma::vec& e, bool usesvd);

#endif
