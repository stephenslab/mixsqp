#include "objective.h"
#include "misc.h"

using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
double compute_objective_helper (const vec& x, const vec& u, const vec& w, 
				 const vec& z, double lambda);

// FUNCTION DEFINITIONS
// --------------------
// Compute the value of the (unmodified) objective at x.
double compute_objective (const mat& L, const vec& w, const vec& x,
			  const vec& z, double lambda, const vec& e) {
  vec u = L*x + e;
  return compute_objective_helper(x,u,w,z,lambda);
}

// Compute the value of the (unmodified) objective at x, using either
// the exact matrix L, or a low-rank SVD approximation of L.
double compute_objective (const mat& L, const mat& U, const mat& V,
			  const vec& w, const vec& x, const vec& z,
			  double lambda, const vec& e, bool usesvd) {
  vec u;
  if (usesvd)
    u = U*(trans(V)*x);
  else
    u = L*x;
  u += e;
  return compute_objective_helper(x,u,w,z,lambda);
}

// Compute the value of the (unmodified) objective at x, after
// precalculating u = L*x + e.
double compute_objective_helper (const vec& x, const vec& u, const vec& w, 
				 const vec& z, double lambda) {
  if (u.min() <= 0)
    stop("Objective is -Inf");
  return lambda*norm2(x)/2 - sum(w % (z + log(u)));
}  
