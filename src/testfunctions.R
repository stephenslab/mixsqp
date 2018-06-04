
REBayes_default = function(L){
  res = KWDual(L, rep(1,dim(L)[2]), rep(1,dim(L)[1])/dim(L)[1])
  x = res$f; x[x < 1e-3] = 0; x = x/sum(x)
  return(x)
}

julia_init()
jDo('using LowRankApprox')
jDo('include("../src/mixSQP.jl")')

mixSQP_julia = function(L){
  r2j(L,'L');
  jDo('x_julia = mixSQP(L,lowrank = "qr")');
  x_julia = j2r('x_julia');
  return (x_julia)
}

mixSQP_rcpp = function(L){
  x0 = rep(1,dim(L)[2])/dim(L)[2];
  mixSQP(L, x0, convtol, ptol, eps, sptol, maxiter, maxqpiter, verbose)$x
}

mixSQP_QP = function(L){
  x0 = rep(1,dim(L)[2])/dim(L)[2];
  qrfact = Matrix::qr(L);
  Q = qr.Q(qrfact)[,1:qrfact$rank]
  R = qr.R(qrfact)[1:qrfact$rank,order(qrfact$pivot)]
  mixSQP_qp(Q, R, x0, convtol, ptol, eps, sptol, maxiter, maxqpiter, verbose)$x
}

mixSQP_R = function(L){
  x0 = rep(1,dim(L)[2])/dim(L)[2];
  mixSQP_r(L, x0, convtol, ptol, eps, sptol, maxiter, maxqpiter, verbose)$x
}


simple_test <- function(n,m, t = 100){
  check_for_objective <- function(values){
    tol <- 1e-6
    eval_f <- function(x){-sum(log(L%*%x + 1e-8))/n}
    max_error <- max(c(eval_f(values[[1]])/eval_f(values[[2]]),
                       eval_f(values[[2]])/eval_f(values[[3]]),
                       eval_f(values[[3]])/eval_f(values[[4]]),
                       eval_f(values[[4]])/eval_f(values[[5]]),
                       eval_f(values[[5]])/eval_f(values[[1]])
    ))
    abs(max_error-1) < tol
  }
  set.seed(2018)
  L = testdata(n,m, type = 2)
  x0 = rep(1,m)/m
  out = microbenchmark("REBayes" = {x_rebayes = REBayes_default(L)},
                       "mixSQP_R" = {x_r = mixSQP_R(L)},
                       "mixSQP_julia" = {x_julia = mixSQP_julia(L)},
                       "mixSQP_Rcpp" = {x_rcpp = mixSQP_rcpp(L)},
                       "mixSQP_QP" = {x_qp = mixSQP_QP(L)}
                       #,check = check_for_objective
                       ,
                       times = t
  )
  return (out)
}
