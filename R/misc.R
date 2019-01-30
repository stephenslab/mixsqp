# Verify that a logical argument is either TRUE or FALSE.
verify.logical.arg <- function (x, arg.name = deparse(substitute(x))) {
  arg.name <- sprintf("\"%s\"",arg.name)
  if (!(is.atomic(x) &
        is.logical(x) &
        length(x) == 1 &
        all(!is.na(x)) &
        all(x == TRUE | x == FALSE)))
    stop(paste("Argument",arg.name,"should be TRUE or FALSE"))
  return(TRUE)
}

# Verify that a non-negative scalar argument is satisfactory.
verify.nonneg.scalar.arg <- function (x, arg.name = deparse(substitute(x))) {
  arg.name <- sprintf("\"%s\"",arg.name)
  if (!(is.atomic(x) &
        is.numeric(x) &
        length(x) == 1 &
        all(!is.na(x)) &
        all(is.finite(x)) &
        all(x >= 0)))
    stop(paste("Argument",arg.name,"should be a non-negative number"))
  return(TRUE)
}

# Verify that a "maxiter" argument---that is, an argument giving the
# maximum number of iterations---is valid. It should be a positive,
# finite, non-missing integer.
verify.maxiter.arg <- function (x, arg.name = deparse(substitute(x))) {
  arg.name <- sprintf("\"%s\"",arg.name)
  if (!(is.atomic(x) &
        is.numeric(x) &
        length(x) == 1 &
        all(!is.na(x)) &
        all(is.finite(x)) &
        all(x > 0) &
        all(round(x) == x)))
    stop(paste("Argument",arg.name,"should be an integer value",
               "greater than zero"))
  return(TRUE)
}

# Verify that the likelihood matrix specifying the optimization
# problem is valid. The likelihood matrix should be a numeric matrix
# with at least two columns, and all the entries should be positive.
# It is assumed that the input argument is named "L".
#
# If the matrix is not valid, an error is reported; otherwise, TRUE is
# returned.
verify.likelihood.matrix <- function (L) {
  msg <- paste("Input argument \"L\" should be a numeric matrix with >= 2",
               "columns, >= 1 rows, all its entries should be non-negative,",
               "finite and not NA, and some entries should be positive")
  if (!is.matrix(L))
    stop(msg)
  else if (!(nrow(L) >= 1 &
             ncol(L) >= 1 &
             is.numeric(L)))
    stop(msg)
  else if (!(all(L >= 0) &
             all(is.finite(L)) &
             !any(is.na(L)) &
             any(L > 0)))
    stop(msg)
  return(TRUE)
}

# Verify that the vector weights specifying the optimization problem
# is valid, then return the normalized weights in double-precision. It
# is assumed that the weights argument is named "w", and that the
# likelihood matrix argument is named "L". The weights should be a
# numeric vector with all non-negative entries, in which the length is
# equal to the number of rows of L. Furthermore, the weights should
# sum to 1; if not, the weights must be normalized to sum to 1.
#
# Input L should be the provided likelihood matrix.
# 
# If the weights are not valid, an error is reported; otherwise, the
# normalized weights (coerced to double-precision) are returned.
verify.weights <- function (L, w) {
  msg <- paste("Input argument \"w\" should be a numeric vector with",
               "non-negative, finite and non-missing entries, and with",
               "one entry per row of L")
  if (!(is.atomic(w) &
        is.numeric(w)))
    stop(msg)
  else if (!(all(w >= 0) &
             all(is.finite(w)) &
             !any(is.na(w)) &
             length(w) == nrow(L)))
    stop(msg)
  storage.mode(w) <- "double"
  return(w/sum(w))
}

# Verify that the estimate of the solution to the optimization problem
# is valid, then normalize this estimate if necessary. Argument x
# should be a numeric vector with non-negative entries, in which the
# length is equal to the number of columns of L.
#
# Input L should be the provided likelihood matrix. It is assumed that
# the argument providing L is named "L".
#
# If x is not valid, an error is reported; otherwise, the normalized
# initial estimate (coerced to double-precision) is returned.
verify.estimate <- function (x, L, arg.name = deparse(substitute(x))) {
  arg.name <- sprintf("\"%s\"",arg.name)
  msg <- paste("Argument",arg.name,"should be a numeric vector with",
               "only non-negative, finite and non-missing entries, with one",
               "entry per column of L")
  if (!(is.atomic(x) &
        is.numeric(x)))
    stop(msg)
  if (!(all(x >= 0) &
        all(is.finite(x)) &
        !any(is.na(x)) &
        length(x) == ncol(L)))
    stop(msg)
  storage.mode(x) <- "double"
  return(x/sum(x))
}

# Generates a vector of n points that are equally spaced on the
# logarithmic scale. Note that x and y should be positive numbers.
logspace <- function (x, y, n)
  2^seq(log2(x),log2(y),length = n)

# Scale each column A[,i] by b[i].
scale.cols <- function (A, b) {
    
  # TO DO: Modify this code to avoid the transpose of A, which could
  # be a large matrix.
  t(t(A) * b)
}
