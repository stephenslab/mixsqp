
REBayes = function(L){
  res = KWDual(L, rep(1,dim(L)[2]), rep(1,dim(L)[1])/dim(L)[1])
  x = res$f; x[x < 1e-3] = 0; x = x/sum(x)
  return(x)
}

julia_init()
jDo('using LowRankApprox')
jDo('include("../mixsqp-paper/code/julia/mixSQP.jl")')
mixSQP_julia = function(L){
  r2j(L,'L')
  jDo('x_julia = mixSQP(L,lowrank = qr, verbose = false)["x"]')
  x_julia = j2r('x_julia')
  return (x_julia)
}


simple_test <- function(n,m){
  check_for_objective <- function(values){
    tol <- 1e-6
    eval_f <- function(x){-sum(log(L%*%x + 1e-8))/n}
    max_error <- max(c(eval_f(values[[1]])/eval_f(values[[2]]),
                       eval_f(values[[2]])/eval_f(values[[3]]),
                       eval_f(values[[1]])/eval_f(values[[1]])
    ))
    abs(max_error-1) < tol
  }
  L = testdata(n,m)
  x0 = rep(1,m)/m
  out = microbenchmark("REBayes" = {x_rebayes = REBayes(L)},
                       "mixSQP_julia" = {x_julia = mixSQP_julia(L)},
                       "mixSQP_Rcpp" = {x_rcpp = mixSQP_rcpp(L)},
                       check = check_for_objective
  )
  return (out)
}