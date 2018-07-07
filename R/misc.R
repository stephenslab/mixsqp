# Generates a vector of n points that are equally spaced on the
# logarithmic scale. Note that x and y should be positive numbers.
logspace <- function (x, y, n)
  2^seq(log2(x),log2(y),length = n)

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

# Verify that the likelihood matrix specifying the optimization
# problem is valid. The likelihood matrix should be a numeric matrix
# with at least two columns, and all the entries should be
# positive.
#
# If the matrix is not valid, an error is reported; otherwise, TRUE is
# returned.
verify.likelihood.matrix <- function (L) {
  if (!is.matrix(L))
    stop("Argument \"L\" should be a matrix; see \"help(matrix)\"")
  if (!(ncol(L) >= 2 & is.numeric(L) & all(L > 0)))
    stop(paste("Input \"L\" should be a numeric matrix with >= 2 columns,",
               "and all its entries should be positive"))
  return(TRUE)
}

# Verify that the vector weights specifying the optimization problem
# is valid, then return the normalized weights. The weights should be
# a numeric vector with all positive entries, in which the length is
# equal to the number of rows of L. Further, the weights should sum to
# 1;  if not, the weights are normalized to sum to 1.
#
# Input L should be the provided likelihood matrix. If the weights are
# not valid, an error is reported; otherwise, the normalized weights
# (coerced to double-precision) are returned.
verify.weights <- function (L, w) {
    
  if (!(is.vector(w) & is.numeric(w)))
    stop(paste("Argument \"w\" should be a numeric vector;",
               "for more information, see \"help(is.vector)\""))
  if (!(length(w) == nrow(L) & all(w > 0)))
    stop(paste("Input vector \"w\" should contain only positive values,",
               "with one entry per row of L"))
  storage.mode(w) <- "double"
  return(w/sum(w))
}
