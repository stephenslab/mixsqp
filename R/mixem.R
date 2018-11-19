mixem <- function (L, w, x0, convtol = 1e-6, maxiter = 1000, verbose = TRUE) {
  x <- x0
  for (iter in 1:maxiter) {

    # Save the current estimate of the mixture weights and the current
    # objective function value.
    f0 <- f
    x0 <- x

    # E STEP
    # Compute the posterior probabilities
    P <- scale.cols(L,x)
    P <- P / (rowSums(P) + eps)

    # M STEP
    # Update the mixture weights.
    x <- colMeans(P)
    
    # CHECK CONVERGENCE
    # Print the status of the algorithm and check the convergence
    # criterion. Convergence is reached when the maximum difference
    # between the mixture weights at two successive iterations is less
    # than the specified tolerance, or when objective increases.
    maxd <- max(abs(w - w0))
    if (verbose)
      cat(sprintf("%4d %0.3e\n",iter,maxd))
    if (maxd < convtol)
      break
  }
  return(x)
}
