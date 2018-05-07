
mixSQP = function(L, x = rep(1,dim(L)[2]), lowrank = "svd",
                  eps = 1e-8, convtol = 1e-10, sptol = 1e-3, ptol = 1e-12, rank = 20,
                  maxiter = 100, maxqpiter = 100){
  # make x sum up to 1
  x = x/sum(x)
  
  # Get the number of rows (n) and columns (m) of L
  n = dim(L)[1]; m = dim(L)[2];
  
  # low-rank approximation
  t_lowrank = 0
  if (lowrank == "qr_julia"){
    t_lowrank = system.time(f <- qr_julia(L, ptol = ptol))[3]
  } else if (lowrank == "svd_julia"){
    t_lowrank = system.time(f <- svd_julia(L, ptol = ptol))[3]
  } else if (lowrank == "svd"){
    t_lowrank = system.time(f <- svd_R(L, rank = rank))[3]
  }
  names(t_lowrank) = NULL
  
  # timer
  t_gradhess = 0
  t_activeset = 0
  t_linesearch = 0
  
  # start loop
  for (i in 1:maxiter){
    
    # timer
    t1 = Sys.time()
    
    # compute objective gradient hessian
    if (lowrank == "svd"){
      D = as.vector( 1 / (f$Q %*% ( t((t(f$u) * f$d[1:rank])) %*% (t(f$v) %*% x)) + eps))
      G = (f$Q %*% f$u %*% (t(f$v) * f$d[1:rank])) * D;
    } else if(lowrank == "qr_julia"){
      D = as.vector( 1 / (f$Q %*% (f$R %*% x) + eps))
      G = (f$Q * D) %*% f$R;
    } else if(lowrank == "svd_julia"){
      D = as.vector( 1 / ( f$u %*% (diag(f$d) %*% (t(f$v) %*% x)) + eps))
      G = ( f$u %*% (t(f$v) * f$d)) * D;
    } else{
      D = as.vector(1/(L %*% x + eps));
      G = L*D;
    }
    g = -colSums(G) / n;
    H = t(G) %*% G / n + eps * diag(m);
    
    # timer
    t2 = Sys.time()
    
    # Check convergence of outer loop
    if(min(g + 1) >= -convtol) break;
    
    # Initialize the solution to the QP subproblem (y).
    if (i > 1){
      ind    = which(x > sptol);
      y      = rep(0,m);
      y[ind] = 1/length(ind);
    } else{
      ind       = c(1,m);
      y         = rep(0,m);
      y[c(1,m)] = 0.5;
    }
    
    for(j in 1:maxqpiter){
      # Define the smaller QP subproblem.
      s   = length(ind);
      H_s = H[ind,ind];
      d   = H %*% y + 2*g + 1;
      d_s = d[ind];
      
      # Solve the smaller problem.
      p      = rep(0,m);
      p_s    = -solve(H_s,d_s);
      p[ind] = p_s;
      
      # If reached at the solution in the current active set
      if (sqrt(sum(p_s^2)) < convtol){
        if (all(d >= -convtol)){
          break; # Check global convergence using KKT
        } else{
          ind_min = which.min(d);          # Find an index with smallest multiplier
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
    
    # timer
    t3 = Sys.time()
    
    # Perform backtracking line search
    for (t in 1:10){
      if (lowrank == "svd"){
        D_new_inv = as.vector(f$Q %*% ( t((t(f$u) * f$d[1:rank])) %*% (t(f$v) %*% y)) + eps)
      } else if(lowrank == "qr_julia"){
        D_new_inv = as.vector(f$Q %*% (f$R %*% y) + eps)
      } else if(lowrank == "svd_julia"){
        D_new_inv = as.vector(f$u %*% (diag(f$d) %*% (t(f$v) %*% y)) + eps)
      } else{
        D_new_inv = as.vector(L %*% y + eps);
      }
      if (sum(log(D)) + sum(log(D_new_inv)) > sum((x-y) * g) / 2) break;
      y = (y-x)/2 + x;
    }
    
    # timer
    t4 = Sys.time()
    
    # Update the solution to the original optimization problem.
    x = y;
    
    # update timer
    t_gradhess = t_gradhess + t2-t1
    t_activeset = t_activeset + t3-t2
    t_linesearch = t_linesearch + t4-t3
  }
  x[x < sptol] = 0; x = x/sum(x);
  return(list(x = x, numiter = i,
              t = c(t_lowrank = t_lowrank, t_gradhess = t_gradhess,
                    t_activeset = t_activeset, t_linesearch = t_linesearch)))
}