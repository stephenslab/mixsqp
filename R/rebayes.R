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
#' missing). For large matrices, it is preferrable that the matrix is
#' stored in double-precision; see \code{\link{storage.mode}}.
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

  # When all the entries of one or more columns are zero, the mixture
  # weights associated with those columns are necessarily zero. Here
  # we handle this situation.
  nonzero.cols <- which(apply(L,2,max) > 0)
  if (length(nonzero.cols) == 1) {
    warning(paste("All columns of \"L\" are zeros except one; this",
                  "corresponds to the trivial solution \"x\" in which",
                  "x[i] = 1 for one column i, and all other entries of",
                  "\"x\" are zero. No optimization algorithm was needed."))
    x               <- rep(0,m)
    x[nonzero.cols] <- 1
    names(x)        <- colnames(L)
    return(list(x      = x,
                value  = mixobj(L,w,x),
                status = NULL))
  } else if (length(nonzero.cols) < m) {
    warning(paste("One or more columns of \"L\" are all zeros; solution",
                  "entries associated with these columns are trivially",
                  "zero"))
    L <- L[,nonzero.cols]
  }
  
  # SOLVE OPTIMIZATION PROBLEM USING MOSEK
  # --------------------------------------
  # Check that the REBayes package is available.
  if(!requireNamespace("REBayes",quietly = TRUE))
    stop("mixKWDual requires package REBayes")
  d   <- rep(1,ncol(L))
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

  # Compute the value of the objective at the estimated solution.
  f <- mixobj(L,w,x)
  
  # If necessary, insert the zero mixture weights associated with the
  # columns of zeros.
  if (length(nonzero.cols) < m) {
    xnz <- x
    x   <- rep(0,m)
    x[nonzero.cols] <- xnz
  }
  
  # Return a list containing the solution to the optimization problem
  # (x), the value of the objective at the solution (value), and the
  # MOSEK convergence status (status).
  names(x) <- colnames(L)
  return(list(x = x,value = f,status = out$status))
}
