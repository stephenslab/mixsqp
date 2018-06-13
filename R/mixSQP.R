#' @title mixSQP
#' @description mixSQP solves a convex optimization problem from nonparametric maximum-likelihood estimation of mixture proportions. It implements a sequential quadratic programming with active-set subproblem solvers. For gigantic data, use low-rank approximation to speed up the computation.
#' When L is a (n) by (m) matrix of nonnegative entries, mixSQP maximizes
#' the following objective function
#' \deqn{f(x) = \sum_j w_j log (\sum_k L_jk x_k)}
#' subject to the (unit) probability simplex constraint
#' \deqn{\sum_k x_k = 1, x_k \ge 0}
#' Without loss of generality \eqn{\sum_j w_j = 1} is required.
#' The problem is originally from nonparametric empirical Bayes mixture MLE problem (NPMLE or NPEB).
#' @param L a matrix of log-likelihoods of mixture components (n by m)
#' @param x0 a initial value for the optimization problem (default rep(1,m)/m).
#' @param w a vector of weight on each data point (default rep(1,n)/n).
#' @param optmethod a programming language used for solving the problem c("Rcpp","R")
#' @param lowrank a type of low-rank approximation c("none","qr","svd")
#' @param lowrankmethod determines what library is used for low-rank approximation c("R_matrix","Julia_lowrankapprox")
#' @param lowranktol a tolerance used for low-rank approximation (default 1e-5): for enough accuracy, set at most 1e-4 for "R_matrix" and 1e-10 for "Julia_lowrankapprox"
#' @param convtol a convergence tolerance used for algorithm's convergence criterion
#' @param sparsetol a tolerance used for determining active indices
#' @param eps a small constant to safeguard from a numerical issue (default 1e-6).
#' @param maxiter a maximum number of outer loop iterations, determining how many qp subproblems will be solved at most.
#' @param maxqpiter a maximum number of inner loop iterations, determining how many active-set subproblems will be solved at most.
#' @param verbose a logical indicating if it shows progress of the algorithm at each iteration
#' @return returns a solution x (in the current version).
#' @examples
#' n = 1e5; m = 2e1;
#' L = testdata(n,m) # create some simulated data
#' x0 = rep(1,m)/m; # initialization
#' w = rep(1,n)/n; # weight
#' optmethod = "Rcpp"; lowrank = "qr"; lowrankmethod = "R_Matrix";
#' mixSQP(L, x0, optmethod, lowrank, lowrankmethod); # using default tolerances
#' @useDynLib mixSQP
#' @importFrom Rcpp evalCpp
#' @export
mixSQP = function(L, x0 = rep(1,dim(L)[2])/dim(L)[2], w = rep(1,dim(L)[1])/dim(L)[1],
                  optmethod = "Rcpp", lowrank = "none",
                  lowrankmethod = "Julia_lowrankapprox", lowranktol = 1e-10, 
                  convtol = 1e-8, sparsetol = 1e-3, eps = 1e-6,
                  maxiter = 50, maxqpiter = 100, verbose = T){
  
  # TO DO : match.arg
  
  if ((lowrankmethod == "Julia_lowrankapprox") & (lowrank == "qr")){
    if (!require(rjulia))
      stop("no rjulia package installed: you must install rjulia package or use another
           low-rank approximaton method (arg: lowrankmethod)")
    
    require("rjulia");
    if (j2r('typeof(Pkg.installed("LowRankApprox")) == Void')){
      cat("installing LowRankApprox package in Julia")
      jDo('Pkg.add("LowRankApprox")')
    }
    if (!(j2r('isdefined(:LowRankApprox)')) ){
      jDo('using LowRankApprox');
    }
  }
  
  # timer
  t1 = Sys.time();
  
  if (lowrank == "qr"){
    if (lowrankmethod == "R_Matrix"){
      qrfact = Matrix::qr(L, tol = lowranktol);
      Q = qr.Q(qrfact)[,1:qrfact$rank]
      R = qr.R(qrfact)[1:qrfact$rank,order(qrfact$pivot)]
      cat()
    } else if(lowrankmethod == "Julia_lowrankapprox"){
      r2j(L,"L"); r2j(lowranktol,"lowranktol");
      jDo("Q,R,P = pqr(L, rtol = lowranktol[1]); R = R[:,sortperm(P)]");
      Q = j2r("Q"); R = j2r("R");
    } else{
      stop("Error : lowrank:", lowrank," does not support ",lowrankmethod," option.")
    }
  }
  
  # timer
  t2 = Sys.time();
  
  if (optmethod == "Rcpp"){
    if (lowrank == "none"){
      out = mixSQP_rcpp_noapprox(L, x0, w, convtol, sparsetol, eps, maxiter, maxqpiter, verbose)
    } else if (lowrank == "qr"){
      out = mixSQP_rcpp_qr(Q, R, x0, w, convtol, sparsetol, eps, maxiter, maxqpiter, verbose)
    } else{
      stop("Error : optmethod:", optmethod," does not support ",lowrank," option.")
    }
  } else if (optmethod == "R"){
    if (lowrank == "none"){
      out = mixSQP_r_noapprox(L, x0, w, convtol, sparsetol, eps, maxiter, maxqpiter, verbose)
    } else if (lowrank == "qr"){
      out = mixSQP_r_qr(Q, R, x0, w, convtol, sparsetol, eps, maxiter, maxqpiter, verbose)
    } else{
      stop("Error : optmethod:", optmethod," does not support ",lowrank," option.")
    }
  } else{
    stop("Error : optmethod:", optmethod," is an invalid option")
  }
  
  # timer
  t3 = Sys.time();
  
  # show timings
  if (verbose){
    if (lowrank == "qr")
      cat("A low-rank approximation using",lowrank,"took",t2-t1,"seconds\n");
    cat("A convex programming took",t3-t2,"seconds\n");
  }
  
  return (list(x = as.vector(out$x_sparse), niter = out$niter,
               status = ifelse(out$niter < maxiter, "OPTIMAL", "SUBOPTIMAL")))
}