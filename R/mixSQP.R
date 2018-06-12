
#' mixSQP solves a convex optimization problem originated from
#' 
#' (Adaptive SHrinkage, see https://github.com/stephens999/ashr)
#' When L is a (n) by (m) matrix of nonnegative entries, mixSQP maximizes
#' the objective function
#' f(x) = \sum_{j=1}^n w_j \log  \sum_{k=1}^m L_{jk} x_{k}
#' subject to the (unit) probability simplex constraint
#' \sum_{k=1}^m x_{k} = 1, x_k >= 0
#' under additional constraint \sum_{j=1}
#' 
#' \sum_{j=1}^n \log \sum_{k=1}^m L_{jk} x_{k} + \sum_{k=1}^m w_{k} \log x_{k}
#' @param L a matrix of log-likelihoods of mixture components
#' @param x0 a initial value for the optimization problem
#' @param optmethod
#' @param outputlevel controls a level of output
#' @return returns a list of 
#' @examples
#' n = 1e4; m = 1e1;
#' L = testdata(n,m) # create some simulated data
#' x0 = rep(1,m)/m;
#' optmethod = "Rcpp"; lowrank = "qr"; lowrankmethod = "R_Matrix";
#' mixSQP(L, x0, optmethod, lowrank, lowrankmethod);
#' @export

mixSQP = function(L, x0 = rep(1,dim(L)[2])/dim(L)[2], optmethod = "Rcpp", lowrank = "none",
                  lowrankmethod = "Julia_lowrankapprox", lowranktol = 1e-10, 
                  convtol = 1e-8, sparsetol = 1e-3, eps = 1e-8,
                  maxiter = 100, maxqpiter = 100, verbose = T){
  
  if ((lowrankmethod == "Julia_lowrankapprox") & (lowrank == "qr")){
    if (!("rjulia" %in% (.packages())) ){
      require("rjulia");
      jDo('using LowRankApprox');
    }
  }
  
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
  
  if (optmethod == "Rcpp"){
    if (lowrank == "none"){
      mixSQP_rcpp_noapprox(L, x0, convtol, sparsetol, eps, maxiter, maxqpiter, verbose)
    } else if (lowrank == "qr"){
      mixSQP_rcpp_qr(Q, R, x0, convtol, sparsetol, eps, maxiter, maxqpiter, verbose)
    } else{
      stop("Error : optmethod:", optmethod," does not support ",lowrank," option.")
    }
  } else if (optmethod == "R"){
    if (lowrank == "none"){
      mixSQP_r_noapprox(L, x0, convtol, sparsetol, eps, maxiter, maxqpiter, verbose)
    } else if (lowrank == "qr"){
      
    } else{
      stop("Error : optmethod:", optmethod," does not support ",lowrank," option.")
    }
  } else{
    stop("Error : optmethod:", optmethod," is an invalid option")
  }
}