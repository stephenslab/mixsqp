mixSQP = function(L, x0 = rep(1,dim(L)[2])/dim(L)[2], optmethod = "Rcpp", lowrank = "none",
                  lowrankmethod = "Julia_lowrankapprox", lowranktol = 1e-10, 
                  convtol = 1e-8, sparsetol = 1e-3, eps = 1e-8,
                  maxiter = 100, maxqpiter = 100, verbose = T){
  
  if (lowrankmethod == "Julia_lowrankapprox"){
    library("rjulia")
    jDo('using LowRankApprox');
  }
  
  if (lowrank == "qr"){
    if (lowrankmethod == "R_matrix"){
      qrfact = Matrix::qr(L, tol = lowranktol);
      Q = qr.Q(qrfact)[,1:qrfact$rank]
      R = qr.R(qrfact)[1:qrfact$rank,order(qrfact$pivot)]
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
      mixSQP_rcpp_qp(Q, P, x0, convtol, sparsetol, eps, maxiter, maxqpiter, verbose)
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