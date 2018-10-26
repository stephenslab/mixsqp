#' @title Maximum-likelihood estimation of mixture proportions using MOSEK
#' 
#' @description This function is simply a wrapper to
#' \code{\link[REBayes]{KWDual}} with a similar interface to
#' \code{mixsqp} for solving the same problem as \code{mixsqp}. See
#' \code{\link{mixsqp}} and \code{\link[REBayes]{KWDual}} for details.
#'
#' @param L Matrix specifying the optimization problem to be solved.
#' In the context of mixture-model fitting, \code{L[j,k]} should be
#' the value of the kth mixture component density at the jth data
#' point. \code{L} should be a numeric matrix with at least two
#' columns, with all entries being non-negative and finite (and not
#' missing). Further, no column should be entirely zeros. For large
#' matrices, it is preferrable that the matrix is stored in
#' double-precision; see \code{\link{storage.mode}}.
#'
#' @param w An optional numeric vector, with one entry for each row of
#' \code{L}, specifying the "weights" associated with the rows of
#' \code{L}. All weights must be finite, non-negative and not
#' missing. Internally, the weights are normalized to sum to 1,
#' which does not change the problem, but does change the value of the
#' objective function reported. By default, all weights are equal.
#' 
#' @param ... Additional arguments passed to \code{\link[REBayes]{KWDual}}.
#'
#' @return A list object with the following elements:
#'
#' \item{x}{The estimated solution to the convex optimization problem.}
#'
#' \item{value}{The value of the objective function at \code{x}.}
#'
#' \item{status}{The return status from MOSEK.}
#'
#' For more information on this output, see
#' \code{\link[REBayes]{KWDual}} and \code{\link[Rmosek]{mosek}}.
#'
#' @seealso \code{\link{mixsqp}}, \code{\link[REBayes]{KWDual}}
#' 
#' @export
#' 
mixkwdual <- function (L, w = rep(1,nrow(L)), ...)  {

  # CHECK & PROCESS INPUTS
  # ----------------------
  # Check input L and, if necessary, coerce the likelihood matrix to
  # be in double precision.
  verify.likelihood.matrix(L)
  if (storage.mode(L) != "double")
    storage.mode(L) <- "double"

  # Get the number of rows (n) and columns (m) of the matrix L.
  n <- nrow(L)
  m <- ncol(L)

  # Check and process the weights.
  w <- verify.weights(L,w)

  # SOLVE OPTIMIZATION PROBLEM USING MOSEK
  # --------------------------------------
  # Check that the REBayes package is available.
  if(!requireNamespace("REBayes",quietly = TRUE))
    stop("mixKWDual requires package REBayes")
  d   <- rep(1,m)
  out <- REBayes::KWDual(L,d,w,...)

  # POST-PROCESSING STEPS
  # ---------------------
  # Retrieve the dual solution (which is the solution we are
  # interested in).
  x <- out$f
  
  # Make sure the solution is (primal) feasible. Sometimes the
  # solution contains small negative values, or is not exactly
  # normalized.
  x[x < 0] <- 0
  x        <- x/sum(x)

  # Return a list containing the solution to the optimization problem
  # (x), the value of the objective at the solution (value), and the
  # MOSEK convergence status (status).
  names(x) <- colnames(L)
  return(list(x = x,value = mixobj(L,w,x),status = out$status))
}
