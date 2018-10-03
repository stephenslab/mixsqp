// This is included to suppress the warnings from solve() when the
// system is singular or close to singular.
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>

// This depends statement is needed to tell R where to find the
// additional header files.
//
// [[Rcpp::depends(RcppArmadillo)]]
//

using namespace Rcpp;

// FUNCTION DECLARATIONS
// ---------------------
double mixobjective (const arma::mat& L, const arma::vec& w,
		     const arma::vec& x, double e, arma::vec& u);
void   computegrad  (const arma::mat& L, const arma::vec& w,
		     const arma::vec& x, double e, arma::vec& g,
		     arma::mat& H, arma::vec& u, arma::mat& Z,
		     const arma::mat& I);

// FUNCTION DEFINITIONS
// --------------------
// SQP algorithm for computing a maximum-likelihood estimate of a
// mixture model. For more information, see the help and comments
// accompanying the mixsqp R function.
// 
// [[Rcpp::export]]
List mixSQP_rcpp (const arma::mat& L, const arma::vec& w, const arma::vec& x0, 
                  double convtolsqp, double convtolactiveset,
		  double zerothreshold, double eps, double delta,
		  int maxitersqp, int maxiteractiveset, bool verbose) {
  
  // Get the number of rows (n) and columns (m) of the conditional
  // likelihood matrix.
  int n = L.n_rows;
  int m = L.n_cols;

  // Print a brief summary of the analysis, if requested.
  if (verbose) {
    Rprintf("Running mix-SQP 0.1-29 on %d x %d matrix\n",n,m);
    Rprintf("convergence tol. (SQP):  %0.1e\n",convtolsqp);
    Rprintf("conv. tol. (active-set): %0.1e\n",convtolactiveset);
    Rprintf("max. iter (SQP):         %d\n",maxitersqp);
    Rprintf("max. iter (active-set):  %d\n",maxiteractiveset);
    Rprintf("zero threshold:          %0.1e\n",zerothreshold);
  }
  
  // PREPARE DATA STRUCTURES
  // -----------------------
  // Initialize storage for the outputs obj, gmin, nnz, nqp and dmax.
  arma::vec obj(maxitersqp);
  arma::vec gmin(maxitersqp);
  arma::vec nnz(maxitersqp);
  arma::vec nqp(maxitersqp);
  arma::vec nls(maxitersqp);
  arma::vec dmax(maxitersqp);
  
  // Initialize the solution.
  arma::vec x = x0;
  
  // Initialize storage for matrices and vectors used in the
  // computations below.
  arma::vec  g(m);    // Vector of length m storing the gradient.
  arma::vec  u(n);    // Vector of length n storing L*x + eps or its log.
  arma::mat  H(m,m);  // m x m matrix storing Hessian.
  arma::mat  Z(n,m);  // n x m matrix Z = D*L, where D = diag(1/(L*x+e)).
  arma::mat  I(m,m);  // m x m diagonal matrix e*I.
  arma::uvec t(m);    // Temporary unsigned int. vector result of length m.
  arma::vec  y(m);    // Vector of length m storing y
  arma::vec  p(m);    // Vector of length m storing y
  arma::vec  b(m);    // Vector of length m storing H*y+2+g+1
  arma::vec  d(m);    // Vector of length m storing absolute
		      // differences between between two solution
		      // estimates.
  
  int    newind;        // New index to be added or deleted.
  double alpha  = 1;    // Define step size
  double status = 1;    // Convergence status.
  
  // This is used in computing the Hessian matrix.
  I  = arma::eye(m,m);
  I *= delta;
  
  // Initialize some loop variables used in the loops below.
  int i = 0; 
  int j = 0;
  
  // Print the column labels for reporting the algorithm's progress.
  if (verbose)
    Rprintf("iter    objective max.diff max(rdual) nnz nqp nls\n");
  
  // Repeat until the convergence criterion is met, or until we reach
  // the maximum number of (outer loop) iterations.
  for (i = 0; i < maxitersqp; i++) {
    
    // COMPUTE GRADIENT AND HESSIAN
    // ----------------------------
    computegrad(L,w,x,eps,g,H,u,Z,I);
    
    // Report on the algorithm's progress. Here we compute: the value
    // of the objective at x (obj); the smallest gradient value
    // (gmin), which is used as a convergence criterion; the number of
    // nonzeros in the solution (nnz); and the number of inner-loop
    // iterations (nqp).
    t      = (x >= zerothreshold);
    obj[i] = mixobjective(L,w,x,eps,u);

    // Should be minimum of the nonzero x's only.
    gmin[i] = 1 + g.min();
    nnz[i]  = sum(t);
    if (verbose) {
      if (i == 0)
        Rprintf("%4d %+0.5e       NA %+0.3e%4d  NA  NA\n",
		i + 1,obj[i],-gmin[i],int(nnz[i]));
      else
        Rprintf("%4d %+0.5e %0.2e %+0.3e%4d %3d %3d\n",i + 1,obj[i],
		dmax[i-1],-gmin[i],int(nnz[i]),int(nqp[i-1]),int(nls[i-1]));
    }
    
    // Check convergence.
    //
    // NOTE: I believe -gmin is the max residual ("rdual" on
    // p. 609 of Boyd & Vandenberghe).
    //
    // NOTE: We should only be checking this condition for the nonzero
    // co-ordinates of x.
    //
    if (gmin[i] >= -convtolsqp) {
      status = 0;
      break;
    }

    // Initialize the solution to the QP subproblem (y).
    arma::uvec ind = find(t);
    y.fill(0);
    y.elem(ind).fill(1/nnz[i]);
    
    // Run active set method to solve the QP subproblem.
    for (j = 0; j < maxiteractiveset; j++) {
      
      // Define the smaller QP subproblem.
      b = H*y + 2*g + 1;
      arma::vec bs = b.elem(ind);
      arma::mat Hs = H.elem(ind,ind);
      
      // Solve the smaller problem.
      p.fill(0.0);
      p.elem(ind) = -solve(Hs,bs);
      
      // Reset step size
      alpha = 1;
      
      // Check convergence.
      if (arma::norm(p,2) <= convtolactiveset) {
        
        // Compute the Lagrange multiplier.
        if (b.min() >= -convtolactiveset)
          break;
        
        // Find an index with smallest multiplier, Add this to the
        // inactive set.
        newind     = b.index_min();
        t[newind]  = 1;
        ind        = find(t);
      } else{
        
        // Define step size.
        arma::uvec act = find(p < 0);
        
        if (!act.is_empty()) {
          
          arma::vec alp = -y.elem(act)/p.elem(act);
          newind        = alp.index_min();
          
          if (alp[newind] < 1) {
            
            // Blocking constraint exists; find and delete it.
            alpha          = alp[newind]; 
            t[act[newind]] = 0;
            ind            = find(t);
          }
        }
      }
      
      // Move to the new "inner loop" iterate (y) along the search
      // direction.
      y += alpha * p;
    }
    nqp[i] = j + 1;
    
    // BACKTRACKING LINE SEARCH
    // ------------------------
    // sum(x) = sum(y) = 1 so replacing g by g+1 in dot product of x-y
    // and g has no effect.
    for (j = 0; j < 9; j++) {
      if (obj[i] + sum(log(L * y + eps) % w) > dot(x-y,g)/(2*n)) 
	break;
      y = (y-x)/2 + x;
    }
    nls[i] = j + 1;
    
    // UPDATE THE SOLUTION
    // -------------------
    d       = abs(x - y);
    dmax[i] = d.max();
    x = y;
  }
  
  return List::create(Named("x") = x,Named("status") = status);
}

// Compute the value of the (unmodified) objective at x, assuming x is
// (primal) feasible; arguments L and w specify the objective, and e
// is an additional positive constant near zero. Input argument u is a
// vector of length n used to store an intermediate result used in the
// calculation of the objective.
double mixobjective (const arma::mat& L, const arma::vec& w,
		     const arma::vec& x, double e, arma::vec& u) {
  u = L*x + e;
  return -sum(w % log(u));
}

// Compute the gradient and Hessian of the (unmodified) objective at
// x, assuming x is (primal) feasible; arguments L and w specify the
// objective, and e is an additional positive constant near
// zero. Inputs g and H store the value of the gradient and Hessian (a
// vector of length m, and an m x m matrix). Intermediate results used
// in these calculations are stored in three variables: u, a vector of
// length n; Z, an n x m matrix; and ZW, another n x m matrix.
void computegrad (const arma::mat& L, const arma::vec& w, const arma::vec& x,
		  double e, arma::vec& g, arma::mat& H, arma::vec& u,
		  arma::mat& Z, const arma::mat& I) {
   
  // Compute the gradient g = -L'*u where u = w./(L*x + e), and "./"
  // denotes element-wise division.
  u = L*x + e;
  u = w / u;
  g = -trans(L) * u;
    
  // Compute the Hessian H = L'*U*W*U*L, where U = diag(u) and 
  // W = diag(w), with vector u defined as above.
  Z = L;
  u = L*x + e;
  u = sqrt(w) / u;
  Z.each_col() %= u;
  H = trans(Z) * Z + I;
}

// Find the step size.
// double linesearch (const arma::mat& px, double stepdec, double minstepsize) {
//   double a = 1;  // The candidate step size.

//   // Repeat until the sufficient decrease condition is satisfied, or
//   // when we reach the maximum number of backtracking iterations.
//   while (a >= minstepsize) {

//     // Compute the value of the (modified) objective at the candidate
//     // iterate.
    
//     // The candidate point does not meet our criteria, so decrease
//     // the step size.
//     a = a * stepdec;
//   }
  
//   return alpha;
// }
//     // Compute the value of the (modified) objective at the.

//     // Compute the directional gradient with respect to the (modified)
//     // objective.
//       if (psinew < psi + tau*eta*alpha*dpsi) {
//         x <- xnew
//         break
//       }

//           for (j = 0; j < 9; j++){
//       if (obj[i] + sum(log(L * y + eps) % w) > dot(x-y, g) / (2 * n) ) 
// 	break;
// 	  }
	  
