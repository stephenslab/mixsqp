#' @rdname mixSQP
#' @name mixSQP
#' @export
mixobjective <- function (L, x, w = rep(1,nrow(L)), e = 0) {

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

  # COMPUTE OBJECTIVE
  # -----------------
  return(mixobj(L,w,x,e))
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

