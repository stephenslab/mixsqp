#' @title Compute objective optimized by mixsqp.
#'
#' @description See \code{\link{mixsqp}} for a full description of the
#' objective function optimized by the mix-SQP algorithm.
#' 
#' @param L Matrix specifying the optimization problem to be solved.
#' In the context of mixture-model fitting, \code{L[i,j]} should be
#' the value of the jth mixture component density at the jth data
#' point. \code{L} should be a numeric matrix with at least two
#' columns, with all entries being non-negative and finite (and not
#' missing). Further, no column should be entirely zeros. For large
#' matrices, it is preferrable that the matrix is stored in
#' double-precision; see \code{\link{storage.mode}}.
#'
#' @param x The point at which the objective is evaluated in
#' \code{mixobjective}; see argument \code{x0} in \code{\link{mixsqp}}
#' for details.
#' 
#' @param w An optional numeric vector, with one entry for each row of
#' \code{L}, specifying the "weights" associated with the rows of
#' \code{L}. All weights must be finite, non-negative and not
#' missing. Internally, the weights are normalized to sum to 1,
#' which does not change the problem, but does change the value of the
#' objective function reported. By default, all weights are equal.
#' 
#' @param log When \code{log = TRUE}, the input matrix \code{L} is
#' interpreted as containing the logarithm of the data matrix. 
#' 
#' @return The value of the objective at \code{x}. If any entry of
#' \code{L \%*\% x} is less than or equal to zero, \code{Inf} is
#' returned.
#'
#' @seealso \code{\link{mixsqp}}
#' 
#' @export
#' 
mixobjective <- function (L, x, w = rep(1,nrow(L)), log = FALSE) {

  # CHECK & PROCESS INPUTS
  # ----------------------
  # Check input L and, if necessary, coerce the likelihood matrix to
  # be in double precision.
  verify.logical.arg(log)
  if (!is.matrix(L))
    stop("Input argument \"L\" should be a numeric matrix")
  if (log) {
    out <- normalize.loglikelihoods(L)
    L <- out$L
    s <- out$s
  } else
    s <- rep(1,nrow(L))
  verify.likelihood.matrix(L)
  if (storage.mode(L) != "double")
    storage.mode(L) <- "double"
  
  # Check and process the estimate of the solution.
  x <- verify.estimate(x,L)
  
  # Check and process the weights; the weights are normalized to sum
  # to 1.
  w <- verify.weights(L,w)
  
  # COMPUTE OBJECTIVE VALUE
  # -----------------------
  return(mixobj(L,w,x,s))
}

# Compute the value of the objective at x; arguments L and w specify
# the objective, s is a log-scaling factor for each row of L, and e is
# an additional constant that can be set to a small, positive number
# (zero by default) to better ensure numerical stability of the
# optimization.
mixobj <- function (L, w, x, s = rep(1,nrow(L)), e = 0) {
 y <- drop(L %*% x) + e
 if (all(y > 0))
   return(-sum(w * (s + log(y))))
 else
   return(Inf)
}
