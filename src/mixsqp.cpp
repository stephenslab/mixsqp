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
		     const arma::vec& x, const arma::vec& e, 
		     arma::vec& u);
void   computegrad  (const arma::mat& L, const arma::vec& w,
		     const arma::vec& x, const arma::vec& e, 
		     arma::vec& g, arma::mat& H, arma::vec& u, 
		     arma::mat& Z, const arma::mat& I);
double activesetqp (const arma::mat& H, const arma::vec& g, arma::vec& y,
		    arma::uvec& t, int maxiteractiveset,
		    double zerothresholdsearchdir, 
		    double convtolactiveset);
void backtrackinglinesearch (double f, const arma::mat& L,
			     const arma::vec& w, const arma::vec& g,
			     const arma::vec& x, const arma::vec& p,
			     const arma::vec& eps, double suffdecr,
			     double stepsizereduce, double minstepsize, 
			     double& nls, double& stepsize, arma::vec& y,
			     arma::vec& u);

// FUNCTION DEFINITIONS
// --------------------
// SQP algorithm for computing a maximum-likelihood estimate of a
// mixture model. For more information, see the help and comments
// accompanying the mixsqp R function and, in particular, see how
// mixsqp_rcpp is called inside the mixsqp function.
// 
// [[Rcpp::export]]
List mixsqp_rcpp (const arma::mat& L, const arma::vec& w, const arma::vec& x0, 
                  double convtolsqp, double convtolactiveset,
		  double zerothresholdsolution, double zerothresholdsearchdir,
		  double suffdecr, double stepsizereduce, double minstepsize,
		  const arma::vec& eps, double delta, int maxitersqp,
		  int maxiteractiveset, bool verbose) {
  
  // Get the number of rows (n) and columns (m) of the conditional
  // likelihood matrix.
  int n = L.n_rows;
  int m = L.n_cols;

  // Print a brief summary of the analysis, if requested.
  if (verbose) {
    Rprintf("Running mix-SQP algorithm 0.1-97 on %d x %d matrix\n",n,m);
    Rprintf("convergence tol. (SQP):     %0.1e\n",convtolsqp);
    Rprintf("conv. tol. (active-set):    %0.1e\n",convtolactiveset);
    Rprintf("zero threshold (solution):  %0.1e\n",zerothresholdsolution);
    Rprintf("zero thresh. (search dir.): %0.1e\n",zerothresholdsearchdir);
    Rprintf("l.s. sufficient decrease:   %0.1e\n",suffdecr);
    Rprintf("step size reduction factor: %0.1e\n",stepsizereduce);
    Rprintf("minimum step size:          %0.1e\n",minstepsize);
    Rprintf("max. iter (SQP):            %d\n",maxitersqp);
    Rprintf("max. iter (active-set):     %d\n",maxiteractiveset);
  }
  
  // PREPARE DATA STRUCTURES
  // -----------------------
  // Initialize storage for the outputs obj, gmin, nnz, nqp and dmax.
  arma::vec obj(maxitersqp);
  arma::vec gmin(maxitersqp);
  arma::vec nnz(maxitersqp);
  arma::vec nqp(maxitersqp);
  arma::vec nls(maxitersqp);
  arma::vec stepsize(maxitersqp);
  arma::vec dmax(maxitersqp);
  obj.zeros();
  gmin.zeros();
  nnz.zeros();
  nqp.fill(-1);
  nls.fill(-1);
  stepsize.fill(-1);
  dmax.fill(-1);

  // Initialize the solution.
  arma::vec x = x0;
  
  // Initialize storage for matrices and vectors used in the
  // computations below.
  arma::vec  g(m);    // Vector of length m storing the gradient.
  arma::vec  p(m);    // Vector of length m storing the search direction.
  arma::vec  u(n);    // Vector of length n storing L*x + eps or its log.
  arma::mat  H(m,m);  // m x m matrix storing Hessian.
  arma::mat  Z(n,m);  // n x m matrix Z = D*L, where D = diag(1/(L*x+e)).
  arma::mat  I(m,m);  // m x m diagonal matrix e*I.
  arma::uvec t(m);    // Temporary unsigned int. vector result of length m.
  arma::vec  y(m);    // Vector of length m storing y
  arma::vec  d(m);    // Vector of length m storing absolute
		      // differences between between two solution
		      // estimates.
  
  double status = 1;  // Convergence status.
  
  // This is used in computing the Hessian matrix.
  I  = arma::eye(m,m);
  I *= delta;
  
  // Initialize some loop variables used in the loops below.
  int i = 0; 
  
  // Print the column labels for reporting the algorithm's progress.
  if (verbose) {
    Rprintf("iter        objective max(rdual) nnz stepsize max.diff nqp nls");
    Rprintf("\n");
  }
  
  // Repeat until the convergence criterion is met, or until we reach
  // the maximum number of (outer loop) iterations.
  for (i = 0; i < maxitersqp; i++) {

    // Compute the value of the objective at x (obj).
    obj[i] = mixobjective(L,w,x,eps,u);

    // COMPUTE GRADIENT AND HESSIAN
    // ----------------------------
    computegrad(L,w,x,eps,g,H,u,Z,I);
    
    // Determine the nonzero co-ordinates in the current estimate of
    // the solution, x. This specifies the "inactive set".
    t = (x >= zerothresholdsolution);

    // Report on the algorithm's progress. Here we compute: the
    // smallest gradient value (gmin), which is used as a convergence
    // criterion; and the number of nonzeros in the solution (nnz).
    // Note that only the dual residuals (gmin's) corresponding to the
    // nonzero co-ordinates are relevant.
    gmin[i] = 1 + g.min();
    nnz[i]  = sum(t);
    if (verbose) {
      if (i == 0)
        Rprintf("%4d %+0.9e %+0.3e%4d       NA       NA  NA  NA\n",
		i + 1,obj[i],-gmin[i],int(nnz[i]));
      else
        Rprintf("%4d %+0.9e %+0.3e%4d %0.2e %0.2e %3d %3d\n",i + 1,obj[i],
		-gmin[i],(int) nnz[i],stepsize[i-1],dmax[i-1],
		(int) nqp[i-1],(int) nls[i-1]);
    }
    
    // CHECK CONVERGENCE
    // -----------------
    // Convergence is reached with the maximum dual residual is
    // small. The negative of "gmin" is also the maximum dual residual
    // (denoted as "rdual" on p. 609 of Boyd & Vandenberghe, "Convex
    // Optimization", 2009). Although "gmin" here includes both zero
    // (active) and non-zero (inactive) co-ordinates, this condition
    // is trivially satisfied for the zero co-ordinates as the
    // gradient must be non-negative for these co-ordinates. (See
    // communication with Youngseok on Slack.)
    if (gmin[i] >= -convtolsqp) {
      status = 0;
      i++;
      break;
    }

    // This is also a good point to check for a user interrupt; if the
    // user requests an interrupt, then an exception is thrown and
    // control is returned to the R console.
    Rcpp::checkUserInterrupt();
    
    // SOLVE QUADRATIC SUBPROBLEM
    // --------------------------
    // Run the active-set solver to obtain a search direction.
    nqp[i] = activesetqp(H,g,y,t,maxiteractiveset,zerothresholdsearchdir,
			 convtolactiveset);
    p = y - x;
    
    // BACKTRACKING LINE SEARCH
    // ------------------------
    backtrackinglinesearch(obj[i],L,w,g,x,p,eps,suffdecr,stepsizereduce,
			   minstepsize,nls[i],stepsize[i],y,u);
    
    // UPDATE THE SOLUTION
    // -------------------
    d       = abs(x - y);
    dmax[i] = d.max();
    x = y;
  }

  // CONSTRUCT OUTPUT
  // ----------------
  return List::create(Named("x")         = x,
		      Named("status")    = status,
		      Named("objective") = obj.head(i),
		      Named("max.rdual") = -gmin.head(i),
		      Named("nnz")       = nnz.head(i),
		      Named("stepsize")  = stepsize.head(i),
		      Named("max.diff")  = dmax.head(i),
		      Named("nqp")       = nqp.head(i),
		      Named("nls")       = nls.head(i));
}

// Compute the value of the (unmodified) objective at x, assuming x is
// (primal) feasible; arguments L and w specify the objective, and e
// is an additional positive constant near zero. Input argument u is a
// vector of length n used to store an intermediate result used in the
// calculation of the objective.
double mixobjective (const arma::mat& L, const arma::vec& w,
		     const arma::vec& x, const arma::vec& e, arma::vec& u) {
  u = L*x + e;
  if (u.min() <= 0)
    Rcpp::stop("Halting because the objective function has a non-finite value (logarithms of numbers less than or equal to zero) at the current estimate of the solution");
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
		  const arma::vec& e, arma::vec& g, arma::mat& H, arma::vec& u,
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

// This implements the active-set method from p. 472 of of Nocedal &
// Wright, Numerical Optimization, 2nd ed, 2006.
double activesetqp (const arma::mat& H, const arma::vec& g, arma::vec& y,
		    arma::uvec& t, int maxiteractiveset,
		    double zerothresholdsearchdir, 
		    double convtolactiveset) {
  int    m     = g.n_elem;
  double nnz   = sum(t);
  double alpha = 1;  // The step size.
  int    newind;     // New index to be added or deleted.
  int    j;
  arma::vec b(m);    // Vector of length m storing H*y + 2*g + 1.
  arma::vec p(m);    // Vector of length m storing the search direction.
  
  // Initialize the solution to the QP subproblem, y.
  arma::uvec i0 = find(1 - t);
  arma::uvec i1 = find(t);
  y.fill(0);
  y.elem(i1).fill(1/nnz);
    
  // Run active set method to solve the QP subproblem.
  for (j = 0; j < maxiteractiveset; j++) {
      
    // Define the smaller QP subproblem.
    b = H*y + 2*g + 1;
    arma::vec bs = b.elem(i1);
    arma::mat Hs = H.elem(i1,i1);
      
    // Solve the smaller problem.
    p.fill(0);
    p.elem(i1) = -solve(Hs,bs);
      
    // Reset the step size.
    alpha = 1;
      
    // First check that the search direction is close to zero
    // (according to the "zerothresholdsearchdir" parameter).
    if ((p.max() <= zerothresholdsearchdir) &
	(-p.min() <= zerothresholdsearchdir) &
        (i0.n_elem > 0)) {
        
      // If all the Lagrange multiplers in the working set (that is,
      // zeroed co-ordinates) are positive, or nearly positive, we
      // have reached a suitable solution.
      if (b(i0).min() >= -convtolactiveset) {
	j++;
        break;
      }

      // Find an index with smallest multiplier, and add this to the
      // inactive set.
      newind     = i0[b(i0).index_min()];
      t[newind]  = 1;
      i0         = find(1 - t);
      i1         = find(t);

    // In this next part, we consider adding a co-ordinate to the
    // working set, but only if there is at least 2 non-zero
    // co-ordinates.
    } else if (sum(t) > 1) {
        
      // Define the step size.
      arma::uvec act = find(p < 0);
      if (!act.is_empty()) {
        arma::vec alp = -y.elem(act)/p.elem(act);
        newind        = alp.index_min();
        if (alp[newind] < 1) {
            
          // Blocking constraint exists; find and delete it.
          alpha          = alp[newind]; 
          t[act[newind]] = 0;
          i1             = find(t);
	  i0             = find(1 - t);
        }
      }
    }
      
    // Move to the new "inner loop" iterate (y) along the search
    // direction.
    y += alpha * p;
  }

  return j;
}

// This implements the backtracking line search algorithm from p. 37
// of Nocedal & Wright, Numerical Optimization, 2nd ed, 2006. The
// return value is the number of line search iterations needed to
// identify a step size satisfying the "sufficient decrease"
// condition.
// 
// Note that sum(x) = sum(y) = 1, so replacing g by g + 1 in dot product
// of p and g has no effect.
void backtrackinglinesearch (double f, const arma::mat& L,
			     const arma::vec& w, const arma::vec& g,
			     const arma::vec& x, const arma::vec& p,
			     const arma::vec& eps, double suffdecr,
			     double stepsizereduce, double minstepsize, 
			     double& nls, double& stepsize, arma::vec& y,
			     arma::vec& u) {
  double fnew;
  stepsize = 1;
  nls      = 0;

  // Iteratively reduce the step size until either (1) we can't reduce
  // any more (because we have hit the minimum step size constraint),
  // or (2) the new candidate solution satisfies the "sufficient
  // decrease" condition.
  while (stepsizereduce * stepsize >= minstepsize) {
    y    = x + stepsize*p;
    fnew = mixobjective(L,w,y,eps,u);
    nls++;

    // Check whether the new candidate solution (y) satisfies the
    // sufficient decrease condition. If so, accept this candidate
    // solution.
    if (fnew <= f + suffdecr*stepsize*dot(p,g))
      break;

    // The new candidate does not satisfy the sufficient decrease
    // condition, so we need to try again with a smaller step size.
    stepsize *= stepsizereduce;
  }
}
