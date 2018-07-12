#' @title mixSQP
#' 
#' @description mixSQP solves a convex optimization problem in
#'   https://arxiv.org/abs/1806.01412.
#'   It implements a sequential quadratic programming with
#'   active-set subproblem solvers. For gigantic data, use low-rank
#'   approximation to speed up the computation.
#' 
#' When L is a (n) by (m) matrix of nonnegative entries, mixSQP
#' maximizes the following objective function \deqn{f(x) = \sum_j w_j
#' log (\sum_k L_jk x_k)} subject to the (unit) probability simplex
#' constraint \deqn{\sum_k x_k = 1, x_k \ge 0} Without loss of
#' generality \eqn{\sum_j w_j = 1} is required.
#' 
#' @param L Matrix specifying the optimization problem to be
#'   solved. It should be a numeric matrix with positive entries, and,
#'   ideally, double-precision.
#'
#' @param w A numeric vector, with one entry for each row of \code{L},
#'   specifying the "weights" associated with the rows of \code{L}. All
#'   weights must be positive. It is assumed the weights sum to 1; if
#'   not, they will automatically be normalized to sum to 1. By default,
#'   all rows of \code{L} are assigned the same weight.
#' 
#' @param x0 An initial estimate of the solution to the optimization
#'   problem. It should contain only non-negative values, and the
#'   entries should sum to 1; if not, the vector is automatically
#'   normalized to sum to 1. By default, \code{x0} is the vector with
#'   all equal entries.
#' 
#' @param convtol A convergence tolerance used for algorithm's
#'   convergence criterion.
#' 
#' @param sparsetol A tolerance used for determining active indices.
#' 
#' @param eps A small constant to safeguard from a numerical issue
#'   (default 1e-6).
#' 
#' @param maxiter A maximum number of outer loop iterations,
#'   determining how many qp subproblems will be solved at most.
#' 
#' @param maxqpiter A maximum number of inner loop iterations,
#'   determining how many active-set subproblems will be solved at most.
#' 
#' @param verbose A logical indicating if it shows progress of the
#'   algorithm at each iteration.
#' 
#' @return Returns a solution x (in the current version).
#' 
#' @examples
#' n  <- 1e5
#' m  <- 20
#' w  <- rep(1,n)/n; # weights
#' L  <- testdata(n,m, mix_type = "mix_n") # Create some simulated data
#' x0 <- rep(1,m)/m # initialization
#' out <- mixSQP(L,x0,w)
#' 
#' @useDynLib mixSQP
#' 
#' @importFrom Rcpp evalCpp
#' 
#' @export
#' 
mixSQP <- function(L, w = rep(1,nrow(L)), x0 = rep(1,ncol(L)), 
                   convtol = 1e-8, sparsetol = 1e-8,
                   eps = .Machine$double.eps,
                   maxiter = 1000, maxqpiter = 100, verbose = TRUE){

  # CHECK INPUTS
  # ------------
  # The likelihood matrix should be a numeric matrix with at least
  # two columns, and all the entries should be positive.
  verify.likelihood.matrix(L)

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
  w <- verify.weights(L,w)

  # The initial estimate of the solution should be a numeric vector
  # with all non-negative entries, in which the length is equal to the
  # number fo columns of L. Further, the entries should sum to 1.
  if (!(is.vector(x0) & !is.list(x0) & is.numeric(x0)))
    stop(paste("Argument \"x0\" should be a numeric vector; for more",
               "information, see \"help(is.vector)\" and",
               "\"help(is.numeric)\""))
  if (!(length(x0) == m & all(x0 >= 0)))
    stop(paste("Input vector \"x0\" should contain only non-negative",
               "values, with one entry per column of L"))
  storage.mode(x0) <- "double"
  x0 <- x0/sum(x0)

  # Input argument verbose should be TRUE or FALSE.
  if (!(is.logical(verbose) & length(verbose) == 1))
    stop("Argument \"verbose\" should be TRUE or FALSE")
  
  # SOLVE OPTIMIZATION PROBLEM USING SQP METHOD
  # -------------------------------------------
  out <- mixSQP_rcpp(L,w,x0,convtol,sparsetol,eps,maxiter,maxqpiter,verbose)

  # POST-PROCESSING STEPS
  # ---------------------
  # Label the elements of the solution (x) by the column labels of the
  # likelihood matrix (L).
  x        <- out$x
  x        <- drop(x)
  names(x) <- colnames(L)

  return(list(x = x,value = mixobjective(L,w,x)))
}
