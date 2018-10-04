#' @rdname mixSQP
#' 
#' @name mixSQP
#'
#' @param ... Additional arguments passed to \code{\link[REBayes]{KWDual}}.
#' 
#' @export
#' 
mixKWDual <- function (L, w = rep(1,nrow(L)), ...)  {

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
