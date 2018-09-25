# Verify that the likelihood matrix specifying the optimization
# problem is valid. The likelihood matrix should be a numeric matrix
# with at least two columns, and all the entries should be positive.
# It is assumed that the input argument is named "L".
#
# If the matrix is not valid, an error is reported; otherwise, TRUE is
# returned.
verify.likelihood.matrix <- function (L) {
  msg <- paste("Input argument \"L\" should be a numeric matrix with >= 2",
               "columns, and all its entries should be positive, finite and",
               "not NA")
  if (!is.matrix(L))
    stop(msg)
  else if (!(ncol(L) >= 2 & is.numeric(L)))
    stop(msg)
  else if (!(all(L > 0) & all(is.finite(L)) & !any(missing(L))))
    stop(msg)
  return(TRUE)
}

# Verify that the vector weights specifying the optimization problem
# is valid, then return the normalized weights in double-precision. It
# is assumed that the weights argument is named "w". The weights
# should be a numeric vector with all positive entries, in which the
# length is equal to the number of rows of L. Furthermore, the weights
# should sum to 1; if not, the weights must be normalized to sum to 1.
#
# Input L should be the provided likelihood matrix.
# 
# If the weights are not valid, an error is reported; otherwise, the
# normalized weights (coerced to double-precision) are returned.
verify.weights <- function (L, w) {
  msg <- paste("Input argument \"w\" should be a numeric vector with",
               "positive, finite and non-missing entries, and with one",
               "entry per row of L")
  if (!(is.atomic(w) & is.numeric(w)))
    stop(msg)
  else if (!(all(w > 0) & all(is.finite(w)) & !any(missing(w)) &
             length(w) == nrow(L)))
    stop(msg)
  storage.mode(w) <- "double"
  return(w/sum(w))
}

# Verify that the initial estimate of the solution to the optimization
# problem is valid, then normalize this estimate if necessary. It is
# assumed that the initial estimate argument is named "x0". Argument
# x0 should be a numeric vector with non-negative entries, in which
# the length is equal to the number of columns of L.
#
# Input L should be the provided likelihood matrix. 
#
# If x0 is not valid, an error is reported; otherwise, the normalized
# initial estimate (coerced to double-precision) is returned.
verify.initial.estimate <- function (x0, L) {
  msg <- paste("Argument \"x0\" should be a numeric vector with only",
               "non-negative, finite and non-missing entries, with one",
               "entry per column of L")
  if (!(is.atomic(x0) & is.numeric(x0)))
    stop(msg)
  if (!(all(x0 >= 0) & all(is.finite(x0)) & !any(missing(x0)) &
        length(x0) == ncol(L)))
    stop(msg)
  storage.mode(x0) <- "double"
  return(x0/sum(x0))
}

# Generates a vector of n points that are equally spaced on the
# logarithmic scale. Note that x and y should be positive numbers.
logspace <- function (x, y, n)
  2^seq(log2(x),log2(y),length = n)
