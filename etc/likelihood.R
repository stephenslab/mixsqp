# This file defines functions for computing marginal and conditional
# likelihoods under different models and mixture priors.

# Compute the n x k conditional likelihood matrix, where n is the
# number of samples and k is the number of mixture components, for the
# case when the likelihood is univariate normal and prior is a mixture
# of univariate normals.
condlikmatrix.norm <- function (x, se, s) {

  # Get the number of samples (n) and the number of mixture components (k).
  n <- length(x)
  k <- length(s)

  # Entry (i,j) of the conditional likelihood matrix is equal to
  # N(0,se[i]^2 + s[j]^2), the normal density with zero mean and
  # variance se[i]^2 + s[j]^2.
  L <- matrix(0,n,k)
  for (j in 1:k)
    L[,j] <- dnorm(x,sd = sqrt(se^2 + s[j]^2))
  return(L)
}
