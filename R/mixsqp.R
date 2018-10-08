# Possible convergence status messages in mixsqp.
mixsqp.status.converged      <- "converged to optimal solution"
mixsqp.status.didnotconverge <- "exceeded maximum number of iterations"

#' @title Solution to "Mixture Optimization" Problem
#'
#' @description \code{mixsqp} and \code{mixkwdual} can be used to
#'   compute maximum-likelihood estimates of mixture proportions in a
#'   (finite) mixture model, or it can be used more generally to solve
#'   the a constrained, convex optimization problem of the form given
#'   below (see "Details"). \code{mixsqp} uses a Sequential Quadatric
#'   Programming (SQP) approach to solve a slightly modified primal
#'   formulation of the convex optimization problem. See "References"
#'   for more details about the SQP algorithm and the motivation behind
#'   it. \code{mixkwdual} uses the MOSEK interior-point (IP) algorithm
#'   to solve a dual formulation of the original problem; see
#'   \code{\link[REBayes]{KWDual}} for details.
#'
#' @details Here is a mathematical description of the constrained,
#'   convex optimization problem solved by \code{mixsqp} and
#'   \code{mixkwdual}. Let \eqn{L} be a matrix with \eqn{n} rows and
#'   \eqn{m} rows containing only non-negative entries, and let \eqn{w =
#'   (w_1, \ldots, w_n)} be a matrix of non-negative "weights" that sum
#'   to 1. \code{mixsqp} computes the value of vector \eqn{x = (x_1, \ldots,
#'   x_m)} minimizing the following objective function, \deqn{f(x) =
#'   -1/n \sum_{j=1}^n w_j log (\sum_{k=1}^m L_{jk} x_k),} subject to
#'   the constraint that \eqn{x} lie within the simplex; that is, all
#'   entries of \eqn{x} are non-negative, and the sum of these entries
#'   is equal to 1. The Expectation Maximization (EM) algorithm can be
#'   used to solve this optimization problem, but it is intolerably slow
#'   in many interesting cases.
#'
#'   \code{mixsqp} is implemented using the Armadillo C++ linear
#'   algebra library, which can automatically take advantage of
#'   multithreaded matrix computations to speed up \code{mixsqp} for
#'   large \code{L} matrices, but only when R has been configured with a
#'   multithreaded BLAS/LAPACK library (e.g., OpenBLAS).
#'
#' @param L Matrix specifying the optimization problem to be solved.
#'   In the context of mixture-model fitting, \code{L[j,k]} should be
#'   the response of the kth mixture component density at the jth data
#'   point. \code{L} should be a numeric matrix with at least two
#'   columns, in which all its entries are positive finite (and not
#'   missing). For large matrices, it is preferrable that the matrix is
#'   stored in double-precision; see \code{\link{storage.mode}}.
#'
#' @param w An optional numeric vector, with one entry for each row of
#'   \code{L}, specifying the "weights" associated with the rows of
#'   \code{L}. All weights must be finite, non-negative and not
#'   missing. The weights should sum to 1; if they are not, they will
#'   automatically be normalized to sum to 1. By default, all weights
#'   are the same.
#' 
#' @param x0 An optional numeric vector providing an initial estimate
#'   of the solution to the optimization problem. It should contain only
#'   finite, non-missing, non-negative values, and the entries should
#'   sum to 1; if it is not, the vector is automatically normalized to
#'   sum to 1. By default, \code{x0} is the vector with all equal values.
#' 
#' @param convtol.sqp A small, non-negative number specifying the
#'   convergence tolerance for SQP algorithm; convergence is reached
#'   when the maximum dual residual in the Karush-Kuhn-Tucker optimality
#'   conditions is less than or equal to \code{convtol.sqp}. Smaller
#'   values will result in more stringent convergence criteria and more
#'   accurate solutions, at the expense of greater computation
#'   time. Note that changes to this tolerance parameter may require
#'   respective changes to \code{convtol.activeset} and/or
#'   \code{zero.threshold.searchdir} to obtain reliable convergence.
#'
#' @param convtol.activeset Small, non-negative number specifying the
#'   convergence tolerance for the active-set step. Smaller values will
#'   result in higher quality search directions for the SQP algorithm
#'   but possibly a greater per-iteration computational cost. Note that
#'   changes to this tolerance parameter can affect how reliably the SQP
#'   convergence criterion is satisfied, as determined by
#'   \code{convtol.sqp}.
#' 
#' @param zero.threshold.solution A small, non-negative number used to
#'   determine the "active set"; that is, it determines which entries of
#'   the solution are exactly zero. Any entries that are less than or
#'   equal to \code{zero.threshold.solution} are considered to be
#'   exactly zero. Larger values of \code{zero.threshold.solution} may
#'   lead to speedups for matrices with many columns, at the (slight)
#'   risk of prematurely zeroing some co-ordinates.
#'
#' @param zero.threshold.searchdir A small, non-negative number used
#'   to determine when the search direction in the active-set step is
#'   considered "small enough". Note that changes to this tolerance
#'   parameter can affect how reliably the SQP convergence criterion is
#'   satisfied, as determined by \code{convtol.sqp}, so choose this
#'   parameter carefully.
#' 
#' @param eps A small, non-negative number added to the terms inside the
#'   logarithms to sidestep computing logarithms of zero. This prevents
#'   numerical problems at the cost of introducing a small inaccuracy in
#'   the solution. Increasing this number may lead to faster convergence
#'   but possibly a less accurate solution.
#'
#' @param delta A small, non-negative number added to the diagonal of the
#'   Hessian to improve numerical stability (and possibly the speed)
#'   when computing the search directions in the active-set step. 
#'  
#' @param maxiter.sqp Maximum number of SQP iterations to run before
#' reporting a convergence failure; that is, the maximum number of
#' quadratic subproblems that will be solved by the active-set method.
#' 
#' @param maxiter.activeset Maximum number of active-set iterations
#'   taken to solve each of the quadratic subproblems.
#' 
#' @param verbose If \code{verbose = TRUE}, the algorithm's progress
#'   and a summary of the optimization settings are printed to the
#'   console. The algorithm's progress is displayed in a table with one
#'   row per SQP (outer loop) iteration, and with the following columns:
#'   "iter", the (outer loop) SQP iteration; "objective", the value of
#'   the (unmodified) objective function at the current estimate of the
#'   solution; "max.diff", the maximum difference in the estimates
#'   between two successive iterations; "max(rdual)", the maximum "dual
#'   residual" in the Karush-Kuhn-Tucker (KKT) conditions, which is used
#'   to monitor convergence (see \code{convtol.sqp}); "nnz", the number
#'   of non-zero co-ordinates in the current estimate, as determined by
#'   \code{zero.threshold.solution}, "nqp", the number of (inner loop)
#'   active-set iterations taken to solve the quadratic subproblem;
#'   "nls", the number of iterations in the backtracking line search.
#' 
#' @return \code{mixobjective} returns the value of the (unmodified)
#' objective at \code{x}.
#'
#' \code{mixkwdual} returns a list object with the following
#' list elements:
#'
#' \item{x}{The estimated solution to the convex optimization problem.}
#'
#' \item{value}{The value of the (unmodified) objective function at
#'   \code{x}.}
#'
#' \item{status}{The return status from MOSEK.}
#'
#' For more information on this output, see
#' \code{\link[REBayes]{KWDual}} and \code{\link[Rmosek]{mosek}}.
#'
#' \code{mixsqp} returns a list object with the following list elements:
#'
#' \item{x}{The estimated solution to the convex optimization problem.}
#'
#' \item{value}{The value of the (unmodified) objective function at
#'   \code{x}.}
#'
#' \item{status}{A character string giving the status of the algorithm
#'   upon termination.}
#'
#' \item{data}{A data frame containing more detailed information about
#'   the algorithm's progress. The data frame has one row per SQP
#'   iteration. For an explanation of the columns, see the description
#'   of the \code{verbose} argument above. Missing values (NA's) in the
#'   last row indicate that these quantities were not computed because
#'   convergence was reached before computing them. Also note that the
#'   association of these quantities is slightly different than the
#'   console output when \code{verbose = TRUE} as the console output
#'   shows some quantities that were computed after the convergence
#'   check in the previous iteration. The last entries of max.diff, nqp
#'   and nls may not have been assigned if the SQP algorithm converged
#'   successfully (as indicated by negative values), in which case we
#'   should more appropriately assign them missing values (NA).}
#'
#' @references
#'   Y. Kim, P. Carbonetto, M. Stephens and M. Anitescu (2018). A fast
#'   algorithm for maximum likelihood estimation of mixture proportions
#'   using sequential quadratic programming. arXiv:1806.01412
#'   \url{https://arxiv.org/abs/1806.01412}.
#'
#' @seealso \code{\link[REBayes]{KWDual}}
#' 
#' @examples
#' set.seed(1)
#' n  <- 1e5
#' m  <- 10
#' w  <- rep(1,n)/n
#' L  <- simulatemixdata(n,m)$L
#' out.mixsqp <- mixsqp(L,w)
#' out.kwdual <- mixkwdual(L,w)
#' print(mixobjective(L,out.mixsqp$x,w),digits = 16)
#' print(mixobjective(L,out.kwdual$x,w),digits = 16)
#' 
#' @useDynLib mixsqp
#' 
#' @importFrom Rcpp evalCpp
#' 
#' @export
#' 
mixsqp <- function (L, w = rep(1,nrow(L)), x0 = rep(1,ncol(L)), 
                    convtol.sqp = 1e-8, convtol.activeset = 1e-10,
                    zero.threshold.solution = 1e-6,
                    zero.threshold.searchdir = 1e-8,
                    eps = .Machine$double.eps, delta = 1e-10,
                    maxiter.sqp = 1000, maxiter.activeset = 100,
                    verbose = TRUE) {

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

  # Input arguments "convtol.sqp", "convtol.activeset",
  # "zero.threshold.solution", "zero.threshold.searchdir", and "eps"
  # should be non-negative scalars.  Additionally,
  # "zero.threshold.solution" should be less than 1/m.
  verify.nonneg.scalar.arg(convtol.sqp)
  verify.nonneg.scalar.arg(convtol.activeset)
  verify.nonneg.scalar.arg(zero.threshold.solution)
  verify.nonneg.scalar.arg(zero.threshold.searchdir)
  verify.nonneg.scalar.arg(eps)
  if (zero.threshold.solution >= 1/m)
    stop(paste("Behavior of algorithm will be unpredictable if",
               "zero.threshold > 1/m, where m = ncol(X)"))
  
  # Input argument "verbose" should be TRUE or FALSE.
  verify.logical.arg(verbose)
  
  # SOLVE OPTIMIZATION PROBLEM USING mix-SQP
  # ----------------------------------------
  out <- mixsqp_rcpp(L,w,x0,convtol.sqp,convtol.activeset,
                     zero.threshold.solution,zero.threshold.searchdir,
                     eps,delta,maxiter.sqp,maxiter.activeset,verbose)
  
  # Get the algorithm convergence status. The convention is that
  # status = 0 means that the algorithm has successfully converged to
  # the optimal solution, and a status = 1 means that the algorithm
  # reached the maximum number of iterations before converging to a
  # solution.
  if (out$status == 0)
    status <- mixsqp.status.converged
  else
    status <- mixsqp.status.didnotconverge
  if (verbose) {
    if (out$status == 0)
      cat("Convergence criteria met---optimal solution found.\n")
    else 
      cat("Failed to converge within iterations limit.\n")
  }
  
  # POST-PROCESS RESULT
  # -------------------
  # The last entries of max.diff, nqp and nls may not have been
  # assigned if the SQP algorithm converged successfully (as indicated
  # by negative values), in which case we should more appropriately
  # assign them missing values (NA).
  out$max.diff[out$max.diff < 0] <- NA
  out$nqp[out$nqp < 0]           <- NA
  out$nls[out$nls < 0]           <- NA

  # Label the elements of the solution (x) by the column labels of the
  # likelihood matrix (L).
  x        <- out$x
  x        <- drop(x)
  names(x) <- colnames(L)

  # CONSTRUCT OUTPUT
  # ----------------
  return(list(x      = x,
              status = status,
              value  = mixobj(L,w,x),
              data   = data.frame(objective = out$objective,
                                  max.rdual = out$max.rdual,
                                  nnz       = out$nnz,
                                  max.diff  = out$max.diff,
                                  nqp       = out$nqp,
                                  nls       = out$nls)))
}
