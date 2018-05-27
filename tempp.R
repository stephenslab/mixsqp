X = matrix(rnorm(5*2),nrow = 5)
a = qr(L)



#microbenchmark("SVD_rcpp" = {dcSVD(X)}, "SVD_R" = {svd(X)})
tic("SVD_R"); a = svd(X); toc();
tic("QR_R"); a = qr(X); toc();
r2j(X,'X')
a = jDo('@time pqrfact(X, rtol = 1e-10)');


X = matrix(rnorm(1000,100),nrow = 100)
microbenchmark("qr1" = {base::qr(X)}, "qr2" = {Matrix::qr(X)})

X = as.matrix(read.table("~/Desktop/sample10000x40.txt")); colnames(L) <- NULL
m <- function(.) as(., "matrix")
qX <- Matrix::qr(  X)
drop0(R. <- qr.R(qX), tol=1e-15)
Q. <- qr.Q(qX)


L = testdata(1e5,50, type = 1)
x = mixSQP_julia(L)
x3 = REBayes(L)
convtol = 1e-8; ptol = 1e-10; eps = 1e-8; sptol = 1e-3; maxiter = 100; maxqpiter = 100; verbose = F;
tic("mixSQP_Rcpp"); x_rcpp <- mixSQP_rcpp(L); toc();
setwd("~/git/mixsqp")
sourceCpp("mixSQP.cpp")
sourceCpp("mixSQP_qp.cpp")
library(tictoc)

L = testdata(2 * 1e5,40, type = 2)
x0 = rep(1,dim(L)[2])/dim(L)[2];
tic(); x = mixSQP_julia(L); toc();
tic(); x3 = REBayes(L); toc();
tic(); x1 = mixSQP_rcpp(L); toc();
tic(); x2 = mixSQP_QP(L); toc();

eval_f = function(x){
  -sum(log(L %*% x + eps))/dim(L)[1]
}

c(eval_f(x),eval_f(x1),eval_f(x2),eval_f(x3)) - min(c(eval_f(x),eval_f(x1),eval_f(x2),eval_f(x3)))
