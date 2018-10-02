#' @title mixSQP
#'
#' @description mixSQP solves a convex optimization problem in
#'   https://arxiv.org/abs/1806.01412.
#'   It implements a sequential quadratic programming with
#'   active-set subproblem solvers.
#' 
#' When L is a (n) by (m) matrix of nonnegative entries, mixSQP
#' maximizes the following objective function \deqn{f(x) = \sum_j w_j
#' \log (\sum_k L_jk x_k)} subject to the (unit) probability simplex
#' constraint \deqn{\sum_k x_k = 1, x_k \ge 0} Without loss of
#' generality \eqn{\sum_j w_j = 1} is required.
#' 
#' @param L Matrix specifying the optimization problem to be solved.
#'   It should be a numeric matrix with at least two columns, with
#'   finite and positive entries, and, ideally, stored in
#'   double-precision.
#'
#' @param w A numeric vector, with one entry for each row of \code{L},
#'   specifying the "weights" associated with the rows of \code{L}. All
#'   weights must be finite and positive. The weights should sum to 1;
#'   if not, they will automatically be normalized to sum to 1. By
#'   default, all weights are set to 1, in which case all rows of
#'   \code{L} are assigned the same weight.
#' 
#' @param x0 An initial estimate of the solution to the optimization
#'   problem. It should contain only finite, non-negative values, and
#'   the entries should sum to 1; if not, the vector is automatically
#'   normalized to sum to 1. By default, \code{x0} is the vector with
#'   all equal entries.
#' 
#' @param convtol.sqp A convergence tolerance used for algorithm's
#'   convergence criterion.
#' 
#' @param zero.threshold A tolerance used for determining active indices.
#' 
#' @param eps A small constant to safeguard from a numerical issue.
#' 
#' @param maxiter.sqp Maximum number of SQP iterations; that is, the
#'   maximum number of quadratic subproblems that will be solved by the
#'   active-set method.
#' 
#' @param maxiter.activeset Maximum number of active-set iterations
#'   taken to solve each of the quadratic subproblems.
#' 
#' @param verbose If \code{verbose = TRUE}, print progress of algorithm
#'   to console.
#' 
#' @return Returns a solution x (in the current version). Also returns
#'   the algorithm convergence status.
#' 
#' @examples
#' set.seed(1)
#' n  <- 1e5
#' m  <- 10
#' w  <- rep(1,n)/n
#' L  <- simulatemixdata(n,m)$L
#' out.mixsqp <- mixSQP(L,w)
#' out.kwdual <- mixKWDual(L,w)
#' print(mixobjective(L,out.mixsqp$x,w),digits = 16)
#' print(mixobjective(L,out.kwdual$x,w),digits = 16)
#' 
#' @useDynLib mixSQP
#' 
#' @importFrom Rcpp evalCpp
#' 
#' @export
#' 
mixSQP <- function(L, w = rep(1,nrow(L)), x0 = rep(1,ncol(L)), 
                   convtol.sqp = 1e-8, zero.threshold = 1e-8,
                   eps = .Machine$double.eps, maxiter.sqp = 1000,
                   maxiter.activeset = 100, verbose = TRUE){

  # CHECK & PROCESS INPUTS
  # ----------------------
  # Check the likelihood matrix and, if necessary, coerce the
  # likelihood matrix to be in double-precision.
  verify.likelihood.matrix(L)
  if (storage.mode(L) != "double")
    storage.mode(L) <- "double"

  # Get the number of rows (n) and columns (m) of the matrix L.
  n <- nrow(L)
  m <- ncol(L)

  # Check and process the weights.
  w <- verify.weights(L,w)

  # Check and process the initial estimate of the solution.
  x0 <- verify.estimate(x0,L)
  
  # Input arguments "maxiter.sqp" and "maxiter.activeset" should be
  # scalars that are integers greater than zero.
  verify.maxiter.arg(maxiter.sqp)
  verify.maxiter.arg(maxiter.activeset)
  maxiter.sqp       <- as.double(maxiter.sqp)
  maxiter.activeset <- as.double(maxiter.activeset)

  # Input arguments "convtol.sqp", "zero.threshold" and "eps" should be
  # non-negative scalars.
  verify.nonneg.scalar.arg(convtol.sqp)
  verify.nonneg.scalar.arg(zero.threshold)
  verify.nonneg.scalar.arg(eps)
  
  # Input argument "verbose" should be TRUE or FALSE.
  verify.logical.arg(verbose)
  
  # SOLVE OPTIMIZATION PROBLEM USING mix-SQP
  # ----------------------------------------
  out <- mixSQP_rcpp(L,w,x0,convtol.sqp,zero.threshold,eps,maxiter.sqp,
                     maxiter.activeset,verbose)

  # Get the algorithm convergence status. Currently, 
  if (out$status == 0) {
    status <- "converged to optimal solution"
    if (verbose)
      cat("Convergence criteria met---optimal solution found.\n")
  } else if (out$status == 1) {
    status <- "exceeded maximum number of iterations"
    if (verbose)
      cat(paste("Failed to converge within iterations limit.\n"))
  }
  
  # POST-PROCESS RESULT
  # -------------------
  # Label the elements of the solution (x) by the column labels of the
  # likelihood matrix (L).
  x        <- out$x
  x        <- drop(x)
  names(x) <- colnames(L)

  return(list(x = x,status = status,value = mixobj(L,w,x)))
}
