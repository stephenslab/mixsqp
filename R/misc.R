# Compute the value of the objective at x; arguments L and w specify
# the objective, and "eps" is an additional constant that can be set
# to a small, positive number (zero by default) to better ensure
# numerical stability of the optimization.
mixobjective <- function (L, w, x, eps = 0) {
 if (any(x < 0))
   return(Inf)
 else
   return(-sum(w * log(drop(L %*% x) + eps)))
}

# Generates a vector of n points that are equally spaced on the
# logarithmic scale. Note that x and y should be positive numbers.
logspace <- function (x, y, n)
  2^seq(log2(x),log2(y),length = n)

