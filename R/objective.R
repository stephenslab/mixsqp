#' @rdname mixSQP
#' @name mixSQP
#' @export
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

