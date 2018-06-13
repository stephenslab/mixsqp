
mixSQP_r_noapprox = function(L, x0 = rep(1,dim(L)[2]), w = rep(1,dim(L)[1]),
                             convtol = 1e-8, sparsetol = 1e-3, eps = 1e-6,
                             maxiter = 50, maxqpiter = 100, verbose = T){
  # make x sum up to 1
  x = x0/sum(x0)
  # make w sum up to 1
  w = w/sum(w);
  
  # Get the number of rows (n) and columns (m) of L
  n = dim(L)[1]; m = dim(L)[2];
  
  # start loop
  for (i in 1:maxiter){
    
    # compute objective gradient hessian
    D = as.vector(1/(L %*% x + eps));
    G = L * D;
    g = -colSums(G * w);
    H = t(G) %*% (G * w) + eps * diag(m);
    
    # Check convergence of outer loop
    if(min(g + 1) >= -convtol) break;
    
    # Initialize the solution to the QP subproblem (y).
    ind    = which(x > sparsetol);
    y      = rep(0,m);
    y[ind] = 1/length(ind);
    
    for(j in 1:maxqpiter){
      # Define the smaller QP subproblem.
      s   = length(ind);
      H_s = H[ind,ind];
      d   = H %*% y + 2*g + 1;
      d_s = d[ind];
      
      # Solve the smaller problem.
      p      = rep(0,m);
      C = chol(H_s);
      p_s    = -backsolve(C,forwardsolve(t(C),d_s))
      p[ind] = p_s;
      
      # If reached at the solution in the current active set
      if (sqrt(sum(p_s^2)) < convtol){
        if (all(d >= -convtol)){
          break; # Check global convergence using KKT
        } else{
          d2 = d; d[ind] = Inf;
          ind_min = which.min(d2);          # Find an index with smallest multiplier
          ind     = sort(c(ind, ind_min)); # Add this to the inactive set
        }
      } else{     
        # Find a feasible step length.
        alpha     = 1;
        alpha0    = -y[ind]/p_s;
        ind_block = which(p_s < 0);
        if (!(length(ind_block) == 0)){
          alpha0  = alpha0[ind_block];
          t       = which.min(alpha0);
          # If there exists a blocking constraint
          if (alpha0[t] < 1){
            # Find it
            alpha = alpha0[t];
            # Update working set
            ind   = ind[-ind_block[t]];
          }
        }
        # Move to the new y along the search direction.
        y = y + alpha * p;
      }
    }
    
    # Perform backtracking line search
    for (t in 1:10){
      D_new_inv = as.vector(L %*% y + eps);
      if (sum(log(D)*w) + sum(log(D_new_inv)*w) > sum((x-y) * g) / (2*n) ) break;
      y = (y-x)/2 + x;
    }
    
    # Update the solution to the original optimization problem.
    x = y;
  }
  x[x < sparsetol] = 0; x = x/sum(x);
  return(list(x_sparse = x, niter = i))
}