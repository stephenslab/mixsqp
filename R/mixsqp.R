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
#' @param convtol A convergence tolerance used for algorithm's
#'   convergence criterion.
#' 
#' @param sparsetol A tolerance used for determining active indices.
#' 
#' @param eps A small constant to safeguard from a numerical issue
#'   (default 1e-6).
#' 
#' @param maxitersqp Maximum number of SQP iterations; that is, the
#'   maximum number of quadratic subproblems that will be solved by the
#'   active-set method.
#' 
#' @param maxiteractiveset Maximum number of active-set iterations
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
                   convtol = 1e-8, sparsetol = 1e-8,
                   eps = .Machine$double.eps,
                   maxitersqp = 1000, maxiteractiveset = 100,
                   verbose = TRUE){

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
  x0 <- verify.estimate(x0,L,"x0")
  
  # Input arguments maxitersqp and maxiteractiveset should both be
  # positive scalars.
  if (!(is.numeric(maxitersqp) & all(maxitersqp >= 1) &
        length(maxitersqp) == 1 & is.numeric(maxiteractiveset) &
        all(maxiteractiveset >= 1) & length(maxiteractiveset) == 1))
    stop(paste("Arguments \"maxitersqp\" and \"maxiteractiveset\" should be",
               "numbers that are 1 or greater"))
  
  # Input argument verbose should be TRUE or FALSE.
  if (!(is.logical(verbose) & length(verbose) == 1))
    stop("Argument \"verbose\" should be TRUE or FALSE")
  
  # SOLVE OPTIMIZATION PROBLEM USING SQP METHOD
  # -------------------------------------------
  out <- mixSQP_rcpp(L,w,x0,convtol,sparsetol,eps,maxitersqp,
                     maxiteractiveset,verbose)

  # Get the algorithm convergence status.
  if (out$status == 0) {
    status <- "converged to optimal solution"
    if (verbose)
      cat("Convergence criteria met; optimal solution found.\n")
  } else if (out$status == 1) {
    status <- "exceeded maximum number of iterations"
    if (verbose)
      cat(paste("Failed to converge within iterations limit.\n"))
  } 
  
  # POST-PROCESSING STEPS
  # ---------------------
  # Label the elements of the solution (x) by the column labels of the
  # likelihood matrix (L).
  x        <- out$x
  x        <- drop(x)
  names(x) <- colnames(L)

  return(list(x = x,status = status,value = mixobj(L,w,x)))
}
