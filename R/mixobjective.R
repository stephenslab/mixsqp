#' Compute the value of the objective at x; arguments L and w specify
#' the objective, and e is an additional constant that can be set to a
#' small, positive number (zero by default) to better ensure numerical
#' stability of the optimization.
#'
#' @export
#' 
mixobjective <- function (L, w = rep(1,nrow(L)), x = rep(1,ncol(L)), e = 0) {
 if (any(x < 0))
   return(Inf)
 else
   return(-sum(w * log(drop(L %*% x) + e)))
}

