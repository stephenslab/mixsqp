#' @title Add title here.
#'
#' @description Add more detailed description here. Explain that this
#'   solves the dual problem formulation using the MOSEK interior-point
#'   solver (see text from manuscript).
#'
#' @param L Matrix specifying the optimization problem to be
#'   solved. It should be a numeric matrix with positive entries, and
#'   ideally double-precision. 
#'
#' @param w A numeric vector, with one entry for each row of \code{L},
#'   specifying the "weights" associated with the rows of \code{L}. All
#'   weights must be positive. It is assumed the weights sum to 1; if
#'   not, they will automatically be normalized to sum to 1. By default,
#'   all rows of \code{L} are assigned the same weight.x
#' 
#' @param ... Additional optimization parameters passed to MOSEK. See
#'   \code{\link[REBayes]{KWDual}} for details.
#'
#' @return \code{mixKWDual} returns a list with two components:
#'
#' \item{x}{The solution to the optimization problem, as provided by 
#'   MOSEK.}
#'
#' \item{value}{The value of the objective at \code{x}.}
#' 
#' \item{status}{The MOSEK convergence status.}
#' 
#' @examples
#' # Add example here.
#'
#' @export
#' 
mixKWDual <- function (L, w = rep(1,nrow(L)), ...)  {

  # CHECK INPUTS
  # ------------
  # The likelihood matrix should be a numeric matrix with at least
  # two columns, and all the entries should be positive.
  if (!is.matrix(L))
    stop("Argument \"L\" should be a matrix; see \"help(matrix)\"")
  if (!(ncol(L) >= 2 & is.numeric(L) & all(L > 0)))
    stop(paste("Input \"L\" should be a numeric matrix with >= 2 columns,",
               "and all its entries should be positive"))

  # If necessary, coerce the likelihood matrix to be in double
  # precision.
  if (storage.mode(L) != "double")
    storage.mode(L) <- "double"

  # Get the number of rows (n) and columns (m) of the matrix L.
  n <- nrow(L)
  m <- ncol(L)

  # The weights should be a numeric vector with all positive entries,
  # in which the length is equal to the number of rows of L. Further,
  # the weights should sum to 1.
  if (!(is.vector(w) & is.numeric(w)))
    stop("Argument \"w\" should be a numeric vector")
  if (!(length(w) == n & all(w > 0)))
    stop(paste("Input vector \"w\" should contain only positive values,",
               "with one entry per row of L"))
  storage.mode(w) <- "double"
  w <- w/sum(w)
  
  # SOLVE OPTIMIZATION PROBLEM USING MOSEK
  # --------------------------------------
  # Check that the REBayes package is available.
  if(!requireNamespace("REBayes",quietly = TRUE))
    stop("mixKWDual requires package REBayes")
  d   <- rep(1,m)
  out <- REBayes::KWDual(L,d,w,...)

  # Retrieve the dual solution (which is the solution we are
  # interested in).
  x <- out$f
  
  # Make sure the solution is (primal) feasible.
  x[x < 0] <- 0
  x        <- x/sum(x)

  # Return a list containing the solution to the optimization problem
  # (x), the value of the objective at the solution (value), and the
  # MOSEK convergence status (status).
  return(list(x = x,value = mixobjective(L,w,x),status = out$status))
}
