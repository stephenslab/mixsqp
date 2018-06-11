#include <RcppArmadillo.h>

// This depends statement is needed to tell R where to find the
// additional header files.
//
// [[Rcpp::depends(RcppArmadillo)]]
//

using namespace Rcpp;

// SQP algorithm for optimizing mixtures. For more information, see
// the help and comments accompanying the "mixsqp" function in R.
// 
// [[Rcpp::export]]
List mixSQP_rcpp_qp   (const arma::mat& Q, const arma::mat& R, const arma::vec& x0,
                      double convtol, double sparsetol, double eps,
                      int maxiter, int maxqpiter,
                      bool verbose) {
  
  // Get the number of rows (n) and columns (k) of the conditional
  // likelihood matrix.
  int n = Q.n_rows;
  int r = Q.n_cols;
  int k = R.n_cols;
  
  // Print a brief summary of the analysis, if requested.
  if (verbose) {
    Rprintf("Running SQP algorithm with the following settings:\n");
    Rprintf("- %d x %d data matrix\n",n,k);
    Rprintf("- convergence tolerance = %0.2e\n",convtol);
    Rprintf("- zero threshold        = %0.2e\n",sparsetol);
  }
  
  // PREPARE DATA STRUCTURES FOR OPTIMIZATION ALGORITHM
  // Initialize storage for the outputs obj, gmin, nnz and nqp.
  arma::vec obj(maxiter);
  arma::vec gmin(maxiter);
  arma::vec nnz(maxiter);
  arma::vec nqp(maxiter);
  
  // Initialize the solution.
  arma::vec x = x0;
  
  // Initialize storage for matrices and vectors used in the
  // computations below.
  arma::vec    g(k);   // Vector of length k storing the gradient.
  arma::vec    u(n);   // Vector of length n storing L*x + eps or its log.
  arma::mat    H(k,k); // k x k matrix storing Hessian.
  arma::mat    Z(n,r); // n x r matrix storing Z
  arma::mat    I(k,k); // k x k diagonal matrix eps*I.
  arma::uvec   t(k);   // Temporary unsigned integer vector result of length k.
  
  arma::vec  y(k);   // Vector of length k storing y
  arma::vec  p(k);   // Vector of length k storing y
  arma::vec  b(k);   // Vector of length k storing H*y+2+g+1
  
  // This is used in computing the Hessian matrix.
  I  = arma::eye(k,k);
  I *= eps;
  
  // Initialize some loop variables used in the loops below.
  int j = 0;
  
  // Print the column labels for reporting the algorithm's progress.
  if (verbose)
    Rprintf("iter       objective -min(g+1) #nnz #qp\n");
  
  // Repeat until we reach the maximum number of outer loop iterations.
  for (int i = 0; i < maxiter; i++) {
    
    // COMPUTE GRADIENT AND HESSIAN
    // Compute u
    u = Q * (R * x) + eps;
    
    // Compute Z
    Z = Q;
    Z.each_col() /= u;
    
    // Compute the gradient g
    g = -R.t() * arma::sum(Z.t(),1)/n;
    
    // Compute the Hessian H
    Z = Q;
    Z.each_col() /= u;
    H = R.t() * (Z.t() * Z) * R/n + I;
    
    // Report on the algorithm's progress. Here we compute: the value
    // of the objective at x (obj); the smallest gradient value
    // (gmin), which is used as a convergence criterion; the number of
    // nonzeros in the solution (nnz); and the number of inner-loop
    // iterations (nqp).
    
    t       = (x > sparsetol);
    obj[i]  = -sum(log(u));
    gmin[i] = 1 + g.min();
    nnz[i]  = sum(t);
    nqp[i]  = j;
    if (verbose)
      Rprintf("%4d %0.8e %+0.2e %4d %3d\n",i,obj[i],-gmin[i],nnz[i],nqp[i]);
    
    // Check convergence.
    
    if (gmin[i] >= 0)
      break;
    
    // Initialize the solution to the QP subproblem (y).
    arma::uvec  ind = find(t);
    y.fill(0);
    y.elem(ind).fill(1/nnz[i]);
    
    // Run active set method to solve the QP subproblem.
    for (j = 0; j < maxqpiter; j++) {
      
      // Define the smaller QP subproblem.
      b = H*y + 2*g + 1;
      arma::vec   bs  = b.elem(ind);
      arma::mat   Hs  = H.elem(ind,ind);
      
      // Solve the smaller problem.
      p.fill(0.0);
      p.elem(ind) = - inv_sympd(Hs) * bs;
      
      int     newind = 0;    // new index to be added or deleted
      double  alpha = 1;     // Define step size
      
      // Check convergence.
      
      if (arma::norm(p,2) < convtol) {
        
        // Compute the Lagrange multiplier.
        if (b.min() >= -convtol)
          break;
        
        // Find an index with smallest multiplier, Add this to the inactive set
        newind     = b.index_min();
        t[newind]  = 1;
        ind        = find(t);
      } else{
        
        // Define step size
        arma::uvec   act = find(p < 0);
        
        if (!act.is_empty()) {
          
          arma::vec  alp = -y.elem(act)/p.elem(act);
          newind         = alp.index_min();
          
          if (alp[newind] < 1) {
            
            // Blocking constraint exists: find and delete it
            alpha          = alp[newind]; 
            t[act[newind]] = 0;
            ind            = find(t);
          }
        }
      }
      
      // Move to the new "inner loop" iterate (y) along the search direction.
      y += alpha * p;
    }
    
    // PERFORM LINE SEARCH
    for (j = 0; j < 10; j++){
      if (obj[i] + sum(log(Q * (R * y) + eps)) > dot(x-y, g)/2 ) break;
      y = (y-x)/2 + x;
    }
    
    // UPDATE THE SOLUTION
    x = y;
  }
  
  arma::vec x_sparse = x;
  x_sparse.elem(find(x < sparsetol)).fill(0);
  
  return List::create(Named("x") = x/sum(x));
}