library("Rcpp")
source("mixSQP.R")
convtol = 1e-8; ptol = 1e-10; eps = 1e-8; sptol = 1e-3;
maxiter = 100; maxqpiter = 100; verbose = TRUE;
sourceCpp("temp.cpp")


# x0 MUST BE "double"
L = as.matrix(read.table("sample5000x20.txt", header = F)); colnames(L) = NULL
n = dim(L)[1]; m = dim(L)[2];
x0 = rep(0,m); x0[c(1,2,m)] = 1; x0 = x0/sum(x0)

system.time(out1 <- mixSQP(L, x0, "none", eps, convtol, sptol, ptol, 0,
                           maxiter, maxqpiter))

system.time(out2 <- mixsqp_rcpp(L, x0, convtol, ptol, eps, sptol, maxiter,
                                maxqpiter, verbose))

out1$x
as.vector(out2$x)

# x0 MUST BE "double"
L = as.matrix(read.table("sample100000x100.txt", header = F));
colnames(L) = NULL
n = dim(L)[1]; m = dim(L)[2];
x0 = rep(0,m); x0[c(1,2,m)] = 1; x0 = x0/sum(x0)

system.time(out1 <- mixSQP(L, x0, "none", eps, convtol, sptol, ptol,
                           0, maxiter, maxqpiter))

system.time(out2 <- mixsqp_rcpp(L, x0, convtol, ptol, eps, sptol, maxiter, maxqpiter, verbose))

out1$x
as.vector(out2$x)

