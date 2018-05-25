X = matrix(rnorm(5*2),nrow = 5)
a = qr(X)



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


L = testdata(1e5,10)
x = mixSQP_julia(L)
convtol = 1e-8; ptol = 1e-10; eps = 1e-7; sptol = 1e-3; maxiter = 100; maxqpiter = 100; verbose = T;
tic("mixSQP_Rcpp"); x_rcpp <- mixSQP_rcpp(L); toc();
