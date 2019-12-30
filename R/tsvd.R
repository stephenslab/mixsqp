# Compute a truncated SVD approximation U*V' to rectangular matrix X,
# such that X is an n x m matrix, U is an n x k matrix, and V is an m
# x k matrix, where k <= min(n,m). The rank, k, is determined by the
# number of singular values surpassing the tolerance, "tol".
#
# Note that this function will only work if min(dim(X)) is 2 or
# greater.
#
# If irlba fails, the return value is NULL.
#
tsvd <- function (X, tol) {
  r <- min(dim(X))
  k <- 2

  # Iteratively increase the number of singular vectors in the SVD
  # until it is accurate enough (based on the tolerance setting,
  # "tol"), or until we hit an upper limit.
  while (TRUE) {
    out <- tryCatch(irlba(X,k,tol = 1,svtol = 0.01*tol),
                    error   = function (e) NULL,
                    warning = function (e) NULL)
    if (is.null(out))
      return(NULL)
    else if (k == r)
      break
    else if (min(out$d) < tol)
      break
    else
      k <- min(2*k,r)
  }

  # Get the truncated SVD.
  i <- which(out$d > tol)
  if (length(i) < 2)
    i <- 1:2
  d <- out$d[i]
  U <- out$u[,i]
  V <- out$v[,i]
  U <- scale.cols(U,sqrt(d))
  V <- scale.cols(V,sqrt(d))
  return(list(U = U,V = V))
}
