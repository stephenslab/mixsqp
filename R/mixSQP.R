#' @title mixSQP
#' @description mixSQP solves a convex optimization problem originated from ASH
#' (Adaptive SHrinkage, see https://github.com/stephens999/ashr)
#' When L is a (n) by (m) matrix of nonnegative entries, mixSQP maximizes
#' the objective function
#' \deqn{f(x) = \sum_j w_j log  \sum_{k=1}^m L_{jk} x_{k}}
#' subject to the (unit) probability simplex constraint
#' \deqn{\sum_k x_k = 1, x_k >= 0}
#' under additional constraint \eqn{\sum_{j=1} w_j = 1}.
#' \deqn{\sum_{j=1}^n log \sum_{k=1}^m L_{jk} x_{k} + \sum_{k=1}^m w_{k} log x_{k}}
#' @param L a matrix of log-likelihoods of mixture components
#' @param x0 a initial value for the optimization problem
#' @param optmethod Describe optmethod here.
#' @param outputlevel controls a level of output
#' @return returns a list of 
#' @examples
#' n = 1e4; m = 1e1;
#' L = testdata(n,m) # create some simulated data
#' x0 = rep(1,m)/m;
#' optmethod = "Rcpp"; lowrank = "qr"; lowrankmethod = "R_Matrix";
#' mixSQP(L, x0, optmethod, lowrank, lowrankmethod);
#' 
#' @useDynLib mixSQP
#' @importFrom Rcpp sourceCpp
#' @export
mixSQP = function(L, x0 = rep(1,dim(L)[2])/dim(L)[2], optmethod = "Rcpp", lowrank = "none",
                  lowrankmethod = "Julia_lowrankapprox", lowranktol = 1e-10, 
                  convtol = 1e-8, sparsetol = 1e-3, eps = 1e-8,
                  maxiter = 50, maxqpiter = 100, verbose = T){
  
  # TO DO : match.arg
  
  if ((lowrankmethod == "Julia_lowrankapprox") & (lowrank == "qr")){
    if (!require(rjulia))
      stop("no rjulia package installed: you must install rjulia package or use another
           low-rank approximaton method (arg: lowrankmethod)")
    
    require("rjulia");
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
      out = mixSQP_rcpp_noapprox(L, x0, convtol, sparsetol, eps, maxiter, maxqpiter, verbose)
    } else if (lowrank == "qr"){
      out = mixSQP_rcpp_qr(Q, R, x0, convtol, sparsetol, eps, maxiter, maxqpiter, verbose)
    } else{
      stop("Error : optmethod:", optmethod," does not support ",lowrank," option.")
    }
  } else if (optmethod == "R"){
    if (lowrank == "none"){
      out = mixSQP_r_noapprox(L, x0, convtol, sparsetol, eps, maxiter, maxqpiter, verbose)
    } else if (lowrank == "qr"){
      out = mixSQP_r_qr(Q, R, x0, convtol, sparsetol, eps, maxiter, maxqpiter, verbose)
    } else{
      stop("Error : optmethod:", optmethod," does not support ",lowrank," option.")
    }
  } else{
    stop("Error : optmethod:", optmethod," is an invalid option")
  }
  
  # timer
  t3 = Sys.time();
  
  # show timings
  if (lowrank == "qr")
    cat("A low-rank approximation using",lowrank,"took",t2-t1,"seconds\n");
  cat("A convex programming took",t3-t2,"seconds\n");
  
  return (out$x)
}