#' @title Compute objective function optimized by mixsqp
#'
#' @details see \code{mixsqp} for details of the objective function
#' 
#' @param L Matrix specifying the optimization problem to be solved.
#'   In the context of mixture-model fitting, \code{L[j,k]} should be
#'   the value of the kth mixture component density at the jth data
#'   point. \code{L} should be a numeric matrix with at least two
#'   columns, with all entries positive and finite (and not
#'   missing). For large matrices, it is preferrable that the matrix is
#'   stored in double-precision; see \code{\link{storage.mode}}.
#'
#' @param x The point at which the objective is evaluated in
#'   \code{mixobjective}; see argument \code{x0} to \code{mixsqp} for details.
#'   
#' @param w An optional numeric vector, with one entry for each row of
#'   \code{L}, specifying the "weights" associated with the rows of
#'   \code{L}. All weights must be finite, non-negative and not
#'   missing. Internally the weights are normalized to sum to one, which does not
#'   change the problem, but does change the value of the objective function reported.  
#'   By default, all weights are equal. 
#' 
#' @export
#' 
mixobjective <- function (L, x, w = rep(1,nrow(L))) {

  # CHECK & PROCESS INPUTS
  # ----------------------
  # Check input L and, if necessary, coerce the likelihood matrix to
  # be in double precision.
  verify.likelihood.matrix(L)
  if (storage.mode(L) != "double")
    storage.mode(L) <- "double"
  
  # Check and process the estimate of the solution.
  x <- verify.estimate(x,L)
  
  # Check and process the weights.
  w <- verify.weights(L,w)
  
  # COMPUTE OBJECTIVE VALUE
  # -----------------------
  return(mixobj(L,w,x))
}

# Compute the value of the objective at x; arguments L and w specify
# the objective, and e is an additional constant that can be set to a
# small, positive number (zero by default) to better ensure numerical
# stability of the optimization.
mixobj <- function (L, w, x, e = 0) {
 if (any(x < 0))
   return(Inf)
 else
   return(-sum(w * log(drop(L %*% x) + e)))
}

