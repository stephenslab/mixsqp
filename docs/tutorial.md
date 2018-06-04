---
title: "mix-SQP"
output: html_document
---

## Test setting

For the test you must be able to install "Rcpp", "rjulia", "REBayes", "Rmosek", "microbenchmark".

This is only required for the following test -- you should be able to run the code itself without any of them.


```r
# set working directory
# setwd("../src")

# install Rcpp
library("Rcpp")
library("RcppArmadillo")

# load mixSQP.cpp
sourceCpp("../src/mixSQP.cpp")
sourceCpp("../src/mixSQP_qp.cpp")
```


```r
# this libraries must be called for the test, but not for the package itself
library("REBayes")
```

```
## Loading required package: Matrix
```

```r
library("rjulia")
library("microbenchmark")
```

```
## Loading required package: microbenchmarkCore
```

```r
# devtools::install_github("armgong/rjulia")
# devtools::install_github("olafmersmann/microbenchmarkCore")
# devtools::install_github("olafmersmann/microbenchmark")
```

You may want to look into the following R files for detailed implementation setting.


```r
# these source codes are only for the test
source("../src/testdata.R")
source("../src/testfunctions.R")
```

```
## Doing 'using LowRankApprox' in a separate thread to force precompilation ...
```

Let's check if the following code is successfully runnable.


```r
# setting
set.seed(2018)
convtol = 1e-8; ptol = 1e-10; eps = 1e-8; sptol = 1e-3; maxiter = 100; maxqpiter = 100; verbose = F;
L = testdata(1e4,20)
print(cbind(x_rebayes = REBayes(L), x_julia = mixSQP_julia(L),
            x_rcpp = mixSQP_rcpp(L), x_rcpp = mixSQP_QP(L)), digits = 3)
```

```
## Error in REBayes(L): could not find function "REBayes"
```

#### Test

For reference, we provide an example to run the test for L of size 1e4 times 20.


```r
# REBayes      : Rmosek version without low-rank approximation
# mixSQP_R     : R version without low-rank approximation
# mixSQP_julia : Julia version without rank revealing QR (RRQR)
# mixSQP_Rcpp  : Rcpp version without low-rank approximation
# mixSQP_QP    : Rcpp version without rank revealing QR (RRQR)
simple_test(1e4,20, t = 10)
```

```
## Error in mixSQP_r(L, x0, convtol, ptol, eps, sptol, maxiter, maxqpiter, : could not find function "mixSQP_r"
```

We also run the test for L of size 1e5 times 100


```r
# REBayes      : Rmosek version without low-rank approximation
# mixSQP_R     : R version without low-rank approximation
# mixSQP_julia : Julia version without rank revealing QR (RRQR)
# mixSQP_Rcpp  : Rcpp version without low-rank approximation
# mixSQP_QP    : Rcpp version without rank revealing QR (RRQR)
simple_test(1e5,100, t = 10)
```

```
## Error in mixSQP_r(L, x0, convtol, ptol, eps, sptol, maxiter, maxqpiter, : could not find function "mixSQP_r"
```

#### Test using other timing packages

For instance, you may want to use another package "tictoc" to count computation times. The following are one example.


```r
# devtools::install_github("jabiru/tictoc")
# library("tictoc")
# n = 1e4; m = 10;
# L = testdata(n,m)
# x0 = rep(1,m)/m
# tic("REBayes"); x_rebayes <- REBayes(L); toc();
# tic("mixSQP_julia"); x_julia <- mixSQP_julia(L); toc();
# tic("mixSQP_Rcpp"); x_rcpp <- mixSQP_rcpp(L); toc();
# tic("mixSQP_QP"); x_rcpp <- mixSQP_QP(L); toc();
```