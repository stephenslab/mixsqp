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
#'   solved. It should be a numeric matrix with positive entries, and
#'   ideally double-precision.
#'
#' @param x0 A initial value for the optimization problem (default rep(1,m)/m).
#' 
#' @param w A vector of weight on each data point (default rep(1,n)/n).
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
mixSQP = function(L, x0 = rep(1/ncol(L),ncol(L)),
                  w = rep(1/nrow(L),nrow(L)),
                  convtol = 1e-8, sparsetol = 1e-3, eps = 1e-6,
                  maxiter = 1000, maxqpiter = 100, verbose = TRUE){

  # (1) CHECK INPUTS
  # ----------------
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
  
  # Input x0 should be a vector of 
    
  out = mixSQP_rcpp(L,x0,w,convtol,sparsetol,eps, maxiter, maxqpiter, verbose)
  return(out)
}
