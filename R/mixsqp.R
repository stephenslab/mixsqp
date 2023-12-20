# Possible convergence status messages in mixsqp.
mixsqp.status.converged      <- "converged to optimal solution"
mixsqp.status.didnotconverge <- "exceeded maximum number of iterations"
mixsqp.status.didnotrun      <- "SQP algorithm was not run"

#' @title Maximum-likelihood estimation of mixture proportions using SQP
#'
#' @description The \code{mixsqp} function uses a Sequential Quadratic
#' Programming (SQP) algorithm to find the maximum likelihood
#' estimates of mixture proportions in a (finite) mixture model. More
#' generally, \code{mixsqp} solves the corresponding constrained,
#' convex optimization problem, which is given below (see
#' \sQuote{Details}). See \sQuote{References} for more details about
#' the SQP algorithm.
#'
#' @details \code{mixsqp} solves the following optimization problem.
#' Let \eqn{L} be a matrix with \eqn{n} rows and \eqn{m} columns
#' containing only non-negative entries, and let \eqn{w = (w_1,
#' \ldots, w_n)} be a vector of non-negative "weights". \code{mixsqp}
#' computes the value of vector \eqn{x = (x_1, \ldots, x_m)}
#' minimizing the following objective function, \deqn{f(x) =
#' -\sum_{j=1}^n w_j \log (\sum_{k=1}^m L_{jk} x_k),} subject to the
#' constraint that \eqn{x} lie within the simplex; that is, the
#' entries of \eqn{x} are non-negative and sum to 1.  Implicitly,
#' there is an additional constraint \eqn{L*x > 0} in order to ensure
#' that the objective has a finite value. In practice, this constraint
#' only needs to be checked for the initial estimate to ensure that it
#' holds for all subsequent iterates.
#' 
#' If all weights are equal, solving this optimization problem
#' corresponds to finding the maximum-likelihood estimate of the
#' mixture proportions \eqn{x} given \eqn{n} independent data points
#' drawn from a mixture model with \eqn{m} components. In this case,
#' \eqn{L_{jk}} is the likelihood for mixture component \eqn{k} and
#' data point \eqn{j}.
#' 
#' The Expectation Maximization (EM) algorithm can be used to solve
#' this optimization problem, but it is intolerably slow in many
#' interesting cases, and mixsqp is much faster. 
#'
#' A special feature of this optimization problem is that the gradient
#' of the objective does not change with re-scaling; for example, if
#' all the entries of matrix \code{L} are multiplied by 100, the
#' gradient does not change. A practical benefit of this property is
#' that the optimization algorithm will behave similarly irrespective
#' of the scale of \code{L}; for example, the same value for the
#' convergence tolerance \code{convtol.sqp} will have the same effect
#' at different scales.
#'
#' A related feature is that the solution to the optimization problem
#' is invariant to rescaling the rows of \code{L}; for example, the
#' solution will remain the same after all the entries in a row of
#' \code{L} are multiplied by 10. A simple normalization scheme
#' divides each row by the largest entry in the row so that all
#' entries of \code{L} are at most 1: \code{L <- L / apply(L,1,max)}
#' Occasionally, it can be helpful to normalize the rows when some of
#' the entries are unusually large or unusually small. This can help
#' to avoid numerical overflow or underflow errors.
#' 
#' The SQP algorithm is implemented using the Armadillo C++ linear
#' algebra library, which can automatically take advantage of
#' multithreaded matrix computations to speed up \code{mixsqp} for
#' large \code{L} matrices, but only when R has been configured with a
#' multithreaded BLAS/LAPACK library (e.g., OpenBLAS).
#'
#' A "debugging mode" is provided to aid in reproducing convergence
#' failures or other issues. When activated, mixsqp will generate an
#' .RData file containing the exact \code{mixsqp} inputs, and will
#' stop execution upon convergence failure. To activate the debugging
#' mode, run \code{options(mixsqp.debug.mode = TRUE)} prior to calling
#' \code{mixsqp}. By default, the output file is \code{mixsqp.RData};
#' the file can be changed by setting the \code{"mixsqp.debug.file"}
#' global option.
#' 
#' The \code{control} argument is a list in which any of the
#' following named components will override the default optimization
#' algorithm settings (as they are defined by
#' \code{mixsqp_control_default}):
#' 
#' \describe{
#'
#' \item{\code{normalize.rows}}{When \code{normalize.rows = TRUE}, the
#' rows of the data matrix \code{L} are automatically scaled so that
#' the largest entry in each row is 1. This is the recommended setting
#' for better stability of the optimization. When \code{log = TRUE},
#' this setting is ignored becase the rows are already normalized.
#' Note that the objective is computed on the original (unnormalized)
#' matrix to make the results easier to interpret.}
#'
#' \item{\code{tol.svd}}{Setting used to determine rank of truncated
#' SVD approximation for L. The rank of the truncated singular value
#' decomposition is determined by the number of singular values
#' surpassing \code{tol.svd}. When \code{tol.svd = 0} or when \code{L}
#' has 4 or fewer columns, all computations are performed using full L
#' matrix.}
#'
#' \item{\code{convtol.sqp}}{A small, non-negative number
#' specifying the convergence tolerance for SQP algorithm; convergence
#' is reached when the maximum dual residual in the Karush-Kuhn-Tucker
#' (KKT) optimality conditions is less than or equal to
#' \code{convtol.sqp}. Smaller values will result in more stringent
#' convergence criteria and more accurate solutions, at the expense of
#' greater computation time. Note that changes to this tolerance
#' parameter may require respective changes to
#' \code{convtol.activeset} and/or \code{zero.threshold.searchdir} to
#' obtain reliable convergence.}
#'
#' \item{\code{convtol.activeset}}{A small, non-negative number
#' specifying the convergence tolerance for the active-set
#' step. Smaller values will result in higher quality search
#' directions for the SQP algorithm but possibly a greater
#' per-iteration computational cost. Note that changes to this
#' tolerance parameter can affect how reliably the SQP convergence
#' criterion is satisfied, as determined by \code{convtol.sqp}.}
#' 
#' \item{\code{zero.threshold.solution}}{A small, non-negative
#' number used to determine the "active set"; that is, it determines
#' which entries of the solution are exactly zero. Any entries that
#' are less than or equal to \code{zero.threshold.solution} are
#' considered to be exactly zero. Larger values of
#' \code{zero.threshold.solution} may lead to speedups for matrices
#' with many columns, at the (slight) risk of prematurely zeroing some
#' co-ordinates.}
#'
#' \item{\code{zero.threshold.searchdir}}{A small, non-negative
#' number used to determine when the search direction in the
#' active-set step is considered "small enough". Note that changes to
#' this tolerance parameter can affect how reliably the SQP
#' convergence criterion is satisfied, as determined by
#' \code{convtol.sqp}, so choose this parameter carefully.}
#'
#' \item{\code{suffdecr.linesearch}}{This parameter determines how
#' stringent the "sufficient decrease" condition is for accepting a
#' step size in the backtracking line search. Larger values will make
#' the condition more stringent. This should be a positive number less
#' than 1.}
#'
#' \item{\code{stepsizereduce}}{The multiplicative factor for
#' decreasing the step size in the backtracking line search. Smaller
#' values will yield a faster backtracking line search at the expense
#' of a less fine-grained search. Should be a positive number less than
#' 1.}
#'
#' \item{\code{minstepsize}}{The smallest step size accepted by the
#' line search step. Should be a number greater than 0 and at most 1.}
#'
#' \item{\code{identity.contrib.increase}}{When the Hessian is not
#' positive definite, a multiple of the identity is added to the
#' Hessian to ensure a unique search direction. The factor for
#' increasing the identity contribution in this modified Hessian is
#' determined by this control parameter.}
#' 
#' \item{\code{eps}}{A small, non-negative number that is added to the
#' terms inside the logarithms to sidestep computing logarithms of
#' zero. This prevents numerical problems at the cost of introducing a
#' small inaccuracy in the solution. Increasing this number may lead
#' to faster convergence but possibly a less accurate solution.}
#'
#' \item{\code{maxiter.sqp}}{Maximum number of SQP iterations to
#' run before reporting a convergence failure; that is, the maximum
#' number of quadratic subproblems that will be solved by the
#' active-set method.}
#' 
#' \item{\code{maxiter.activeset}}{Maximum number of active-set
#' iterations taken to solve each of the quadratic subproblems. If
#' \code{NULL}, the maximum number of active-set iterations is set to
#' \code{min(20,1 + ncol(L))}.}
#'
#' \item{\code{numiter.em}}{Number of expectation maximization (EM)
#' updates to perform prior to running mix-SQP. Although EM can often
#' be slow to converge, this "pre-fitting" step can help to obtain a
#' good initial estimate for mix-SQP at a small cost.}
#' 
#' \item{\code{verbose}}{If \code{verbose = TRUE}, the algorithm's
#' progress and a summary of the optimization settings are printed to
#' the console. The algorithm's progress is displayed in a table with
#' one row per SQP (outer loop) iteration, and with the following
#' columns: "iter", the (outer loop) SQP iteration; "objective", the
#' value of the objective function (see \eqn{f(x)}) at the current
#' estimate of the solution, \eqn{x}; "max(rdual)", the maximum "dual
#' residual" in the Karush-Kuhn-Tucker (KKT) conditions, which is used
#' to monitor convergence (see \code{convtol.sqp}); "nnz", the number
#' of non-zero co-ordinates in the current estimate, as determined by
#' \code{zero.threshold.solution}; "max.diff", the maximum difference
#' in the estimates between two successive iterations; "nqp", the
#' number of (inner loop) active-set iterations taken to solve the
#' quadratic subproblem; "nls", the number of iterations in the
#' backtracking line search.}
#' }
#'
#' @param L Matrix specifying the optimization problem to be solved.
#' In the context of mixture-model fitting, \code{L[j,k]} should be
#' the value of the kth mixture component density at the jth data
#' point. \code{L} should be a numeric matrix with at least two
#' columns, with all entries being non-negative and finite (and not
#' missing). In some cases, it is easier or more natural to compute
#' \code{log(L)}; for example, it is often easier to compute the
#' log-likelihood rather than the likelihood. Setting \code{log = TRUE}
#' will tell \code{mixsqp} to interpret this input as the logarithm of
#' the data matrix. Note that, for large matrices, it is preferrable
#' that the matrix is stored in double-precision; see
#' \code{\link{storage.mode}}.
#'
#' @param w An optional numeric vector, with one entry for each row of
#' \code{L}, specifying the "weights" associated with the rows of
#' \code{L}. All weights must be finite, non-negative and not
#' missing. Internally, the weights are normalized to sum to 1,
#' which does not change the problem, but does change the value of the
#' objective function reported. By default, all weights are equal.
#' 
#' @param x0 An optional numeric vector providing an initial estimate
#' of the solution to the optimization problem. It should contain only
#' finite, non-missing, non-negative values, and all entries of
#' \code{L \%*\% x0} must be greater than zero (to ensure that the
#' objective evaluates to a finite value at \code{x0}). The vector
#' will be normalized to sum to 1. By default, \code{x0} is the vector
#' with all equal values.
#'
#' @param log When \code{log = TRUE}, the input matrix \code{L} is
#' interpreted as containing the logarithm of the data matrix. 
#' 
#' @param control A list of parameters controlling the behaviour of
#' the optimization algorithm. See \sQuote{Details}.
#' 
#' @return A list object with the following elements:
#'
#' \item{x}{If the SQP algorithm converges, this is the solution to
#' the convex optimization problem. If the algorithm fails to
#' converge, it is the best estimate of the solution achieved by the
#' algorithm. Note that if the SQP algorithm terminates before
#' reaching the solution, \code{x} may not satisfy the equality
#' constraint; that is, the entries of \code{x} may not sum to 1.}
#'
#' \item{value}{The value of the objective function, \eqn{f(x)}, at
#' \code{x}.}
#'
#' \item{grad}{The gradient of the objective function at \code{x}.}
#'
#' \item{hessian}{The Hessian of the objective function at
#' \code{x}. The truncated SVD approximation of L is used to compute
#' the Hessian when it is also used for mix-SQP.}
#' 
#' \item{status}{A character string describing the status of the
#' algorithm upon termination.}
#'
#' \item{progress}{A data frame containing more detailed information
#' about the algorithm's progress. The data frame has one row per SQP
#' iteration. For an explanation of the columns, see the description
#' of the \code{verbose} control parameter in \sQuote{Details}. Missing
#' values (\code{NA}'s) in the last row indicate that these quantities were
#' not computed because convergence was reached before computing
#' them. Also note that the storage of these quantities in the
#' \code{progress} data frame is slightly different than in the console
#' output (when \code{verbose = TRUE}) as the console output shows some
#' quantities that were computed after the convergence check in the
#' previous iteration.}
#'
#' @references
#' 
#' Y. Kim, P. Carbonetto, M. Stephens and M. Anitescu (2020). A fast
#' algorithm for maximum likelihood estimation of mixture proportions
#' using sequential quadratic programming. \emph{Journal of
#' Computational and Graphical Statistics} \bold{29},
#' 261-273. \doi{10.1080/10618600.2019.1689985}
#'
#' @seealso \code{\link{mixobjective}}
#' 
#' @examples
#' set.seed(1)
#' n <- 1e5
#' m <- 10
#' w <- rep(1,n)/n
#' L <- simulatemixdata(n,m)$L
#' out.mixsqp <- mixsqp(L,w)
#' f <- mixobjective(L,out.mixsqp$x,w)
#' print(f,digits = 16)
#'
#' @useDynLib mixsqp
#'
#' @importFrom utils modifyList
#' @importFrom utils sessionInfo
#' @importFrom irlba irlba
#' @importFrom Rcpp evalCpp
#' 
#' @export
#' 
mixsqp <- function (L, w = rep(1,nrow(L)), x0 = rep(1,ncol(L)),
                    log = FALSE, control = list()) {

  # SAVE INPUTS (debug mode only)
  # -----------------------------
  if (getOption("mixsqp.debug.mode")) {
    out.file <- getOption("mixsqp.debug.file")
    sinfo    <- sessionInfo()
    message("mixsqp debugging mode is turned on; writing mixsqp inputs and ",
            "sessionInfo to ",out.file)
    save(list = c("L","w","x0","log","control","sinfo"),file = out.file)
  }
    
  # CHECK & PROCESS INPUTS
  # ----------------------
  # Check and process the likelihood matrix and, if necessary, coerce
  # the likelihood matrix to be in double-precision.
  verify.logical.arg(log)
  if (!is.matrix(L))
    stop("Input argument \"L\" should be a numeric matrix")
  if (log)
    L <- normalize.loglikelihoods(L)
  verify.likelihood.matrix(L)
  if (storage.mode(L) != "double")
    storage.mode(L) <- "double"

  # Get the number of rows (n) and columns (m) of the matrix L.
  n      <- nrow(L)
  m      <- ncol(L)
  coords <- colnames(L)
      
  # Check and process the weights.
  w <- verify.weights(L,w)

  # Check and process the initial estimate of the solution. To ensure
  # that the algorithm reaches the solution, we also require that L*x0
  # > 0 hold, which we can check by asking whether the value of the
  # objective is finite.
  x0 <- verify.estimate(x0,L)
  if (is.infinite(mixobj(L,w,x0)))
    stop(paste("Input \"x0\" is not a valid initial estimate; all entries",
               "of the matrix-vector product L %*% x0 should be positive"))

  # Get the optimization algorithm settings.
  if (!is.list(control))
    stop("Argument \"control\" should be a list")
  control0 <- mixsqp_control_default()
  if (any(!is.element(names(control),names(control0))))
    stop("Argument \"control\" contains unknown parameter names")
  control <- modifyList(control0,control,keep.null = TRUE)
  normalize.rows            <- control$normalize.rows
  tol.svd                   <- control$tol.svd
  convtol.sqp               <- control$convtol.sqp
  convtol.activeset         <- control$convtol.activeset
  zero.threshold.solution   <- control$zero.threshold.solution
  zero.threshold.searchdir  <- control$zero.threshold.searchdir
  suffdecr.linesearch       <- control$suffdecr.linesearch
  stepsizereduce            <- control$stepsizereduce
  minstepsize               <- control$minstepsize
  identity.contrib.increase <- control$identity.contrib.increase
  eps                       <- control$eps
  maxiter.sqp               <- control$maxiter.sqp
  maxiter.activeset         <- control$maxiter.activeset
  numiter.em                <- control$numiter.em
  verbose                   <- control$verbose

  # If the maximum number of active-set iterations is set to NULL, set
  # it to be equal to 1 + ncol(L), or 100, whichever is smaller.
  if (is.null(maxiter.activeset))
    maxiter.activeset <- min(20,m + 1)
  
  # Input arguments "maxiter.sqp" and "maxiter.activeset" should be
  # scalars that are integers greater than zero.
  verify.maxiter.arg(maxiter.activeset)
  verify.maxiter.arg(maxiter.sqp)
  maxiter.sqp        <- as.integer(maxiter.sqp)
  maxiter.activeset  <- as.integer(maxiter.activeset)

  # Input arguments "convtol.sqp", "convtol.activeset",
  # "zero.threshold.solution", "zero.threshold.searchdir",
  # "numiter.em", "identity.contrib.increase", and "eps" should be
  # non-negative scalars. Additionally, "zero.threshold.solution"
  # should be less than 1/m. Also, post a warning if eps is within
  # range of the largest value in one of the rows of the matrix L.
  verify.nonneg.scalar.arg(tol.svd)
  verify.nonneg.scalar.arg(convtol.sqp)
  verify.nonneg.scalar.arg(convtol.activeset)
  verify.nonneg.scalar.arg(numiter.em)
  verify.nonneg.scalar.arg(zero.threshold.solution)
  verify.nonneg.scalar.arg(zero.threshold.searchdir)
  verify.nonneg.scalar.arg(suffdecr.linesearch)
  verify.nonneg.scalar.arg(stepsizereduce)
  verify.nonneg.scalar.arg(minstepsize)
  verify.nonneg.scalar.arg(identity.contrib.increase)
  verify.nonneg.scalar.arg(eps)
  if (!(0 < stepsizereduce & stepsizereduce < 1 &
        identity.contrib.increase > 0  &
        suffdecr.linesearch > 0 &
        minstepsize > 0))
    stop(paste("Control parameter \"stepsizereduce\" must be greater than",
               "0 and less than 1, and \"suffdecr.linesearch\",",
               "\"identity.contrib.increase\" and \"minstepsize\" must be",
               "positive"))
  if (zero.threshold.solution >= 1/m)
    stop(paste("Behavior of algorithm will be unpredictable if",
               "zero.threshold > 1/m, where m = ncol(X)"))
  
  # Input arguments "normalize.rows" and "verbose" should be TRUE or
  # FALSE.
  verify.logical.arg(normalize.rows)
  verify.logical.arg(verbose)
  
  # When all the entries of one or more columns are zero, the mixture
  # weights associated with those columns are necessarily zero. Here
  # we handle this situation.
  nonzero.cols <- which(apply(L,2,max) > 0)
  if (m == 1 | length(nonzero.cols) == 1) {
    warning(paste("Only one column of \"L\" has positive entries; this",
                  "corresponds to the trivial solution \"x\" in which",
                  "x[i] = 1 for one column i, and all other entries of",
                  "\"x\" are zero. No optimization algorithm was needed."))
    x               <- rep(0,m)
    x[nonzero.cols] <- 1
    names(x)        <- coords
    return(list(x        = x,
                status   = mixsqp.status.didnotrun,
                value    = mixobj(L,w,x),
                progress = NULL))
  } else if (length(nonzero.cols) < m) {
    warning(paste("One or more columns of \"L\" are all zeros; solution",
                  "entries associated with these columns are trivially",
                  "zero"))
    L  <- L[,nonzero.cols]
    x0 <- x0[nonzero.cols]
    x0 <- x0/sum(x0)
    m0 <- m
    m  <- length(nonzero.cols)
  } else
    m0 <- m

  # If requested, normalize the rows of L. Note that the rows will
  # already be normalized when log = TRUE.
  if (normalize.rows & !log) {
    out <- normalize.rows(L)
    L   <- out$A
    z   <- log(out$z)
    rm(out)
  } else
    z <- rep(0,n)
  
  # Print a brief summary of the analysis, if requested.
  if (verbose) {
    cat(sprintf("Running mix-SQP algorithm 0.3-54 on %d x %d matrix\n",n,m))
    cat(sprintf("convergence tol. (SQP):     %0.1e\n",convtol.sqp))
    cat(sprintf("conv. tol. (active-set):    %0.1e\n",convtol.activeset))
    cat(sprintf("zero threshold (solution):  %0.1e\n",zero.threshold.solution))
    cat(sprintf("zero thresh. (search dir.): %0.1e\n",
                zero.threshold.searchdir))
    cat(sprintf("l.s. sufficient decrease:   %0.1e\n",suffdecr.linesearch))
    cat(sprintf("step size reduction factor: %0.1e\n",stepsizereduce))
    cat(sprintf("minimum step size:          %0.1e\n",minstepsize))
    cat(sprintf("max. iter (SQP):            %d\n",maxiter.sqp))
    cat(sprintf("max. iter (active-set):     %d\n",maxiter.activeset))
    cat(sprintf("number of EM iterations:    %d\n",numiter.em))
  }

  # COMPUTE TRUNCATED SVD
  # ---------------------
  U <- matrix(0,n,1)
  V <- matrix(0,m,1)
  use.svd <- FALSE
  if (tol.svd > 0 & m > 4) {
    if (verbose)
      cat(sprintf("Computing SVD of %d x %d matrix.\n",n,m))
    t1  <- proc.time()
    out <- tsvd(L,tol.svd)
    t2  <- proc.time()
    if (is.null(out)) {
      if (verbose)
        cat("Matrix is not low-rank; falling back to full matrix.\n")
    } else {
      if (verbose) {
        cat(sprintf("SVD computation took %0.2f seconds.\n",
                    t2["elapsed"] - t1["elapsed"]))
        cat(sprintf("Rank of matrix is estimated to be %d.\n",ncol(out$U)))
      }
      if (ncol(out$U) < m) {

        # Only use the SVD of L if it might be worthwhile to do so.
        U           <- out$U
        V           <- out$V
        rownames(U) <- rownames(L)
        rownames(V) <- colnames(L)
        L           <- tcrossprod(U,V)
        use.svd     <- TRUE
      }
      rm(out)
    }
  }

  # INITIALIZE SOLUTION
  # -------------------
  x <- x0

  # Adjust the numerical safeguard to accommodate negative entries in
  # the SVD reconstruction of L, or L itself (if we ever allow it).
  if (use.svd)
    eps <- eps - min(0,min(tcrossprod(U,V)))
  else
    eps <- eps - min(0,min(L))
  eps <- rep(eps,n)
  
  # RUN A FEW ITERATIONS OF EM
  # --------------------------
  # Print the column labels for reporting the algorithm's progress.
  if (verbose)
    cat("iter        objective max(rdual) nnz stepsize max.diff nqp nls\n")
  t1 <- proc.time()
  if (numiter.em > 0) {
    out <- run.mixem.updates(L,w,x,z,numiter.em,eps,zero.threshold.solution,
                             verbose)
    x <- out$x
    progress.em <- out$progress
    rm(out)
  } else
    progress.em <- NULL
  
  # SOLVE OPTIMIZATION PROBLEM USING mix-SQP
  # ----------------------------------------
  runem <- TRUE
  out <- mixsqp_rcpp(L,U,V,w,z,x,use.svd,runem,convtol.sqp,convtol.activeset,
                     zero.threshold.solution,zero.threshold.searchdir,
                     suffdecr.linesearch,stepsizereduce,minstepsize,
                     identity.contrib.increase,eps,maxiter.sqp,
                     maxiter.activeset,verbose)
  t2 <- proc.time()
  if (verbose)
    cat(sprintf("Optimization took %0.2f seconds.\n",
                t2["elapsed"] - t1["elapsed"]))
 
  # Make sure solution sums to 1.
  x <- drop(out$x)
  x <- x/sum(x)
  
  # Get the algorithm convergence status. The convention is that
  # status = 0 means that the algorithm has successfully converged to
  # the optimal solution, and a status = 1 means that the algorithm
  # reached the maximum number of iterations before converging to a
  # solution.
  if (out$status == 0) {
    status <- mixsqp.status.converged
    if (verbose)
      cat("Convergence criteria met---optimal solution found.\n")
  } else {
    status <- mixsqp.status.didnotconverge
    msg <- paste(strwrap(paste("Failed to converge within iterations limit.",
      "If \"maxiter.sqp\" is small, consider increasing it. Otherwise,",
      "convergence failure is typically a numerical issue remedied by",
      "increasing \"eps\" slightly, at the cost of slightly less accurate",
      "solution; see help(mixsqp). An issue report may also be submitted",
      "to https://github.com/stephenslab/mixsqp/issues, accompanied by an",
      ".rds or .RData file containing the mixsqp inputs. If these inputs",
      "are not accessible, an .RData file containing the inputs can be",
      "generated by setting options(mixsqp.debug.mode = TRUE) before",
      "running mixsqp.")),collapse = "\n")  
    if (getOption("mixsqp.debug.mode"))
      stop(msg)
    else
      warning(msg)
  }

  # POST-PROCESS RESULT
  # -------------------
  # The last entries of stepsize, max.diff, nqp and nls may not have been
  # assigned if the SQP algorithm converged successfully (as indicated
  # by negative values), in which case we should more appropriately
  # assign them missing values (NA).
  out$stepsize[out$stepsize < 0] <- NA
  out$max.diff[out$max.diff < 0] <- NA
  out$nqp[out$nqp < 0]           <- NA
  out$nls[out$nls < 0]           <- NA

  # Compute the objective, gradient and Hessian at the estimated
  # solution.
  e <- control$eps
  f <- mixobj(L,w,x,z)
  u <- drop(L %*% x + e)
  g <- drop((-w/u) %*% L)
  if (use.svd)
    H <- V %*% crossprod(sqrt(w)/u * U) %*% t(V)
  else
    H <- crossprod(sqrt(w)/u * L)

  # If necessary, insert the zero mixture weights associated with the
  # columns of zeros.
  if (m < m0) {
    xnz <- x
    x   <- rep(0,m0)
    x[nonzero.cols] <- xnz
  }
  
  # Label the elements of the solution (x) by the column labels of the
  # likelihood matrix (L).
  names(x)    <- coords
  names(g)    <- colnames(L)
  rownames(H) <- colnames(L)
  colnames(H) <- colnames(L)
  
  # CONSTRUCT OUTPUT
  # ----------------
  return(list(x        = x,
              status   = status,
              value    = f,
              grad     = g,
              hessian  = H,
              progress = rbind(progress.em,
                               data.frame(objective = out$objective,
                                          max.rdual = out$max.rdual,
                                          nnz       = out$nnz,
                                          stepsize  = out$stepsize,
                                          max.diff  = out$max.diff,
                                          nqp       = out$nqp,
                                          nls       = out$nls))))
}

#' @rdname mixsqp
#'
#' @export
#' 
mixsqp_control_default <- function()
  list(normalize.rows            = TRUE,
       tol.svd                   = 1e-6,
       convtol.sqp               = 1e-8,
       convtol.activeset         = 1e-10,
       zero.threshold.solution   = 1e-8,
       zero.threshold.searchdir  = 1e-14,
       suffdecr.linesearch       = 0.01,
       stepsizereduce            = 0.75,
       minstepsize               = 1e-8,
       identity.contrib.increase = 10,
       eps                       = 1e-8,
       maxiter.sqp               = 1000,
       maxiter.activeset         = NULL,
       numiter.em                = 10,
       verbose                   = TRUE)

# This function is used within mixsqp to run several EM updates.
run.mixem.updates <- function (L, w, x, z, numiter, eps, zero.threshold,
                               verbose) {
  out      <- mixem_rcpp(L,w,z,x,eps,numiter,zero.threshold,verbose)
  progress <- data.frame(objective = drop(out$objective),
                         max.rdual = rep(as.numeric(NA),numiter),
                         nnz       = drop(out$nnz),
                         stepsize  = rep(1,numiter),
                         max.diff  = drop(out$max.diff),
                         nqp       = rep(as.numeric(NA),numiter),
                         nls       = rep(as.numeric(NA),numiter))
  return(list(x = drop(out$x),progress = progress))
}
