# TO DO: Add roxygen2 docs here.
mixKWDual <- function (L, ...)  {

  # CHECK INPUTS
  # ------------
  # The likelihood matrix should be a numeric matrix with at least
  # two columns, and all the entries should be positive.
  if (!is.matrix(L))
    stop("Argument \"L\" should be a matrix; see \"help(matrix)\"")
  if (!(ncol(L) >= 2 & all(dat$L > 0)))
    stop(paste("Input matrix \"L\" should have at least columns,",
               "and all its entries should be positive"))

  # Get the number of rows (n) and columns of the matrix L.
  n <- nrow(L)
  m <- ncol(L)
  
  # SOLVE OPTIMIZATION PROBLEM
  # --------------------------
  # Check that the REBayes package is available.
  if(!requireNamespace("REBayes",quietly = TRUE))
    stop("mixKWDual requires package REBayes")
  out <- KWDual(L,rep(1,m),rep(1,n)/n,...)

  # Retrieve the dual solution, which is the solution we are
  # interested in.
  x <- out$f
  
  # Make sure the solution is (primal) feasible.
  x[x < 0] <- 0
  x        <- x/sum(x)

  # Return a list containing the solution to the optimization problem
  # (x) as well as the MOSEK status (status).
  return(list(x = x,status = out$status))
}
