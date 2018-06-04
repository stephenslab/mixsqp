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
source("../src/mixSQP.R")
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
print(cbind(x_rebayes = REBayes_default(L), x_julia = mixSQP_julia(L),
            x_rcpp = mixSQP_rcpp(L), x_rcpp = mixSQP_QP(L)), digits = 3)
```

```
##       x_rebayes  x_julia              
##  [1,]    0.4205 0.420200 0.4205 0.4205
##  [2,]    0.0000 0.000000 0.0000 0.0000
##  [3,]    0.0000 0.000000 0.0000 0.0000
##  [4,]    0.0000 0.000000 0.0000 0.0000
##  [5,]    0.0000 0.000000 0.0000 0.0000
##  [6,]    0.0000 0.000000 0.0000 0.0000
##  [7,]    0.0000 0.000000 0.0000 0.0000
##  [8,]    0.0000 0.000000 0.0000 0.0000
##  [9,]    0.0000 0.000000 0.0000 0.0000
## [10,]    0.0234 0.023409 0.0234 0.0234
## [11,]    0.1713 0.171144 0.1713 0.1713
## [12,]    0.0946 0.094485 0.0946 0.0946
## [13,]    0.2760 0.275840 0.2760 0.2760
## [14,]    0.0121 0.012119 0.0121 0.0121
## [15,]    0.0000 0.000000 0.0000 0.0000
## [16,]    0.0021 0.002098 0.0021 0.0021
## [17,]    0.0000 0.000000 0.0000 0.0000
## [18,]    0.0000 0.000705 0.0000 0.0000
## [19,]    0.0000 0.000000 0.0000 0.0000
## [20,]    0.0000 0.000000 0.0000 0.0000
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
## Unit: milliseconds
##          expr       min        lq      mean    median        uq       max
##       REBayes 232.01567 239.11499 275.55338 249.32965 298.90703 382.12410
##      mixSQP_R  36.42912  41.27385 108.96439 107.91049 174.74801 186.32001
##  mixSQP_julia  26.35281  26.46996  29.10277  28.31316  28.88127  39.62036
##   mixSQP_Rcpp  18.06907  18.10933  19.55340  19.50708  20.20663  23.31385
##     mixSQP_QP  31.84729  33.18994  35.39026  35.25221  36.92327  39.13756
##  neval cld
##     10   c
##     10  b 
##     10 a  
##     10 a  
##     10 a
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
## Unit: milliseconds
##          expr        min         lq       mean     median         uq
##       REBayes 17575.7412 17735.9506 18099.5115 18001.5194 18498.1036
##      mixSQP_R  3174.6504  3264.7944  3348.7197  3365.1095  3424.3652
##  mixSQP_julia   575.1085   599.7174   643.7579   618.6874   659.6526
##   mixSQP_Rcpp  1524.3203  1585.1098  1599.5020  1610.3246  1620.4441
##     mixSQP_QP  3196.1014  3230.1659  3328.3592  3332.7900  3358.9611
##         max neval  cld
##  18893.1629    10    d
##   3488.5559    10   c 
##    846.7327    10 a   
##   1658.4344    10  b  
##   3614.6060    10   c
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
