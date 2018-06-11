mixSQP = function(L, x0 = rep(1,dim(L)[2]), optmethod = "Rcpp", lowrank = "none",
                  lowranktol = 1e-10, 
                  convtol = 1e-8, sparsetol = 1e-3, eps = 1e-8,
                  maxiter = 100, maxqpiter = 100, verbose = T){
  
  if (optmethod == "Rcpp"){
    if (lowrank == "none"){
      mixSQP_rcpp_noapprox(L, x0, convtol, sparsetol, eps, maxiter, maxqpiter, verbose)
    } else if (lowrank == "qp"){
      mixSQP_rcpp_qp(Q, P, x0, convtol, sparsetol, eps, maxiter, maxqpiter, verbose)
    } else{
      stop("Error : optmethod:", optmethod," does not support ",lowrank," option.")
    }
  } else if (optmethod == "R")
    if (lowrank == "none"){
      mixSQP_r_noapprox(L, x0, convtol, sparsetol, eps, maxiter, maxqpiter, verbose)
    } else if (lowrank == "qp"){
      
    } else{
      stop("Error : optmethod:", optmethod," does not support ",lowrank," option.")
    }
}