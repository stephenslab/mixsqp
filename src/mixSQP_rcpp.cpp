#include <RcppArmadillo.h>

// This depends statement is needed to tell R where to find the
// additional header files.
//
// [[Rcpp::depends(RcppArmadillo)]]
//

using namespace Rcpp;

// SQP algorithm for optimizing mixtures. For more information, see
// the help and comments accompanying the mixSQP function in R.
// 
// [[Rcpp::export]]
List mixSQP_rcpp (const arma::mat& L, const arma::vec& w, const arma::vec& x0, 
                  double convtol, double sparsetol, double eps, int maxiter, 
		  int maxqpiter, bool verbose) {
  
  // Get the number of rows (n) and columns (m) of the conditional
  // likelihood matrix.
  int n = L.n_rows;
  int m = L.n_cols;

  // Print a brief summary of the analysis, if requested.
  if (verbose) {
    Rprintf("Running mix-SQP algorithm with the following settings:\n");
    Rprintf("- %d x %d data matrix\n",n,m);
    Rprintf("- convergence tolerance      = %0.1e\n",convtol);
    Rprintf("- zero threshold             = %0.1e\n",sparsetol);
    Rprintf("- outer loop iteration limit = %d\n",maxiter);
    Rprintf("- inner loop iteration limit = %d\n",maxqpiter);
  }
  
  // PREPARE DATA STRUCTURES
  // -----------------------
  // Initialize storage for the outputs obj, gmin, nnz and nqp.
  arma::vec  obj(maxiter);
  arma::vec  gmin(maxiter);
  arma::vec  nnz(maxiter);
  arma::uvec nqp(maxiter);
  arma::uvec nls(maxiter);
  
  // Initialize the solution.
  arma::vec x = x0;
  
  // Initialize storage for matrices and vectors used in the
  // computations below.
  arma::vec  g(m);    // Vector of length m storing the gradient.
  arma::vec  u(n);    // Vector of length n storing L*x + eps or its log.
  arma::mat  H(m,m);  // m x m matrix storing Hessian.
  arma::mat  Z(n,m);  // n x m matrix storing  Z = diag(1/(L*x + eps))*L.
  arma::mat  Zw(n,m); // n x m matrix storing  Z * w
  arma::mat  I(m,m);  // m x m diagonal matrix eps*I.
  arma::uvec t(m);    // Temporary unsigned int. vector result of length m.
  arma::vec  y(m);    // Vector of length m storing y
  arma::vec  p(m);    // Vector of length m storing y
  arma::vec  b(m);    // Vector of length m storing H*y+2+g+1

  int    newind = 0;    // new index to be added or deleted
  double alpha = 1;     // Define step size
  
  // This is used in computing the Hessian matrix.
  I  = arma::eye(m,m);
  I *= eps;
  
  // Initialize some loop variables used in the loops below.
  int i = 0; 
  int j = 0;
  
  // Print the column labels for reporting the algorithm's progress.
  if (verbose)
    Rprintf("iter  objective   -min(g+1) #nnz #nqp #nls\n");
  
  // Repeat until the convergence criterion is met, or until we reach
  // the maximum number of (outer loop) iterations.
  for (i = 0; i < maxiter; i++) {
    
    // COMPUTE GRADIENT AND HESSIAN
    
    // Compute Z = diag(1/(L*x + eps)) * L.
    u = L * x + eps;
    Z = L; 
    Z.each_col() /= u;
    Zw = Z; 
    Zw.each_col() %= w;
    
    // Compute the gradient g = -Z'*1/n.
    g = -arma::sum(Zw.t(),1);
    
    // Compute the Hessian H = Z'*Z/n + eps*I.
    H = Z.t()*Zw + I;
    
    // Report on the algorithm's progress. Here we compute: the value
    // of the objective at x (obj); the smallest gradient value
    // (gmin), which is used as a convergence criterion; the number of
    // nonzeros in the solution (nnz); and the number of inner-loop
    // iterations (nqp).
    t       = (x > sparsetol);
    obj[i]  = -sum(log(u) % w);
    gmin[i] = 1 + g.min();
    nnz[i]  = sum(t);
    if (verbose) {
      if (i == 0)
        Rprintf("%4d  %0.5e %+0.2e %4d \n",i+1,obj[i],-gmin[i],int(nnz[i]));
      else
        Rprintf("%4d  %0.5e %+0.2e %4d %4u %4u\n",i+1,obj[i],-gmin[i],int(nnz[i]),nqp[i-1],nls[i-1]);
    }
    
    // Check convergence.
    if (gmin[i] >= 0)
      break;
    
    // Initialize the solution to the QP subproblem (y).
    arma::uvec ind = find(t);
    y.fill(0);
    y.elem(ind).fill(1/nnz[i]);
    
    // Run active set method to solve the QP subproblem.
    for (j = 0; j < maxqpiter-1; j++) {
      
      // Define the smaller QP subproblem.
      b = H*y + 2*g + 1;
      arma::vec bs  = b.elem(ind);
      arma::mat Hs  = H.elem(ind,ind);
      
      // Solve the smaller problem.
      p.fill(0.0);
      p.elem(ind) = -inv_sympd(Hs) * bs;
      
      // Reset step size
      alpha = 1;
      
      // Check convergence.
      if (arma::norm(p,2) < convtol) {
        
        // Compute the Lagrange multiplier.
        if (b.min() >= -convtol)
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
            
            // Blocking constraint exists: find and delete it
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
    for (j = 0; j < 9; j++){
      if (obj[i] + sum(log(L * y + eps) % w) > dot(x-y, g) / (2 * n) ) 
	break;
      y = (y-x)/2 + x;
    }
    nls[i] = j + 1;
    
    // UPDATE THE SOLUTION
    // -------------------
    x = y;
  }
  
  return List::create(Named("x") = x/sum(x), Named("niter") = i+1);
}
