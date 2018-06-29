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
#' @param ... Additional optimization parameters passed to MOSEK. See
#'   \code{\link[REBayes]{KWDual}} for details.
#'
#' @return \code{mixKWDual} returns a list with two components:
#'
#' \item{x}{The solution to the optimization problems provided by 
#'   MOSEK.}
#'
#' \item{status}{The MOSEK return status.}
#' 
#' @examples
#' # Add example here.
#'
#' @export
#' 
mixKWDual <- function (L, ...)  {

  # CHECK INPUTS
  # ------------
  # The likelihood matrix should be a numeric matrix with at least
  # two columns, and all the entries should be positive.
  if (!is.matrix(L))
    stop("Argument \"L\" should be a matrix; see \"help(matrix)\"")
  if (!(ncol(L) >= 2 & is.numeric(L) & all(dat$L > 0)))
    stop(paste("Input \"L\" should be a numeric matrix with >= 2 columns,",
               "and all its entries should be positive"))

  # If necessary, coerce the likelihood matrix to be in double
  # precision.
  if (storage.mode(L) != "double")
    storage.mode(L) <- "double"

  # Get the number of rows (n) and columns (m) of the matrix L.
  n <- nrow(L)
  m <- ncol(L)
  
  # SOLVE OPTIMIZATION PROBLEM USING MOSEK
  # --------------------------------------
  # Check that the REBayes package is available.
  if(!requireNamespace("REBayes",quietly = TRUE))
    stop("mixKWDual requires package REBayes")
  out <- REBayes::KWDual(L,rep(1,m),rep(1,n)/n,...)

  # Retrieve the dual solution, which is the solution we are
  # interested in.
  x <- out$f
  
  # Make sure the solution is (primal) feasible.
  x[x < 0] <- 0
  x        <- x/sum(x)

  # Return a list containing the solution to the optimization problem
  # (x) as well as the MOSEK status (status).
  return(list(x = x,status = out$status))
}
