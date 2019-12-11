# This is a simple R front-end to the Householder rank-revealing QR
# factorization with column-pivoting as implemented in the Eigen
# library.
#
# X is the matrix to be factorized, and "tol" is the threshold for
# determining which pivots are zero; the rank of the QR decomposition
# is determined by the number of nonzero pivots.
#
# The output is a list containing the factors: the n x k matrix Q, and
# the k x m matrix R, were k is the estimated numerical rank of X.
rrqr <- function (X, tol = 1e-8)
  rrqr_rcpp(X,tol)
