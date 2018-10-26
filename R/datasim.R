#' @title Create likelihood matrix from simulated data set
#' 
#' @description Simulate a data set, then compute the conditional
#' likelihood matrix under a univariate normal likelihood and a
#' mixture-of-normals prior. This models a simple nonparametric
#' Empirical Bayes method applied to simulated data.
#' 
#' @param n Positive integer specifying the number of samples to
#' generate and, consequently, the number of rows of the likelihood
#' matrix L.
#' 
#' @param m Integer 2 or greater specifying the number of mixture
#' components.
#'
#' @param simtype The type of data to simulate. If \code{simtype =
#' "n"}, simulate \code{n} random numbers from a mixture of three
#' univariate normals with mean zero and standard deviation 1, 3 and
#' 6. If \code{simtype = "nt"}, simulate from a mixture of three
#' univariate normals (with zero mean and standard deviations 1, 3 and
#' 5), and a t-distribution with 2 degrees of freedom.
#'
#' @param normalize.rows If \code{normalize.rows = TRUE}, normalize
#' the rows of the likelihood matrix so that the largest entry in each
#' row is 1. The maximum-likelihood estimate of the mixture weights
#' should be invariant to this normalization, and can improve the
#' numerical stability of the optimization.
#' 
#' @return \code{simulatemixdata} returns a list with three list
#' elements:
#'
#' \item{x}{The vector of simulated random numbers (it has length n).}
#'
#' \item{s}{The standard deviations of the mixture components in the
#' mixture-of-normals prior. The rules for selecting the standard
#' deviations are based on the \code{autoselect.mixsd} function from
#' the \code{ashr} package.}
#'
#' \item{L}{The n x m conditional likelihood matrix, in which
#' individual entries (i,j) of the likelihood matrix are given by the
#' normal density function with mean zero and variance \code{1 +
#' s[j]^2}.}
#' 
#' @examples
#' 
#' # Generate the likelihood matrix for a data set with 1,000 samples
#' # and a nonparametric Empirical Bayes model with 20 mixture
#' # components.
#' dat <- simulatemixdata(1000,20)
#'
#' @importFrom stats rnorm
#' @importFrom stats dnorm
#' @importFrom stats rt
#' 
#' @export
#' 
simulatemixdata <- function (n, m, simtype = c("n","nt"), 
                             normalize.rows = TRUE) {

  # CHECK INPUTS
  # ------------
  # Input argument n should be at least 1, and m should be at least 2.
  if (!(is.numeric(n) & n >= 1 & is.finite(n) & !missing(n) &
        round(n) == n & length(n) == 1))
    stop("Argument \"n\" should be a finite, positive integer")
  if (!(is.numeric(m) & m >= 2 & is.finite(m) & !missing(m) &
        round(m) == m & length(m) == 1))
    stop("Argument \"m\" should be a positive integer greater than 1")

  # Get the choice of data to simulate.
  if (!is.character(simtype))
    stop("Argument \"simtype\" should be a character vector")
  simtype <- match.arg(simtype)

  # Input argument normalize.rows should be TRUE or FALSE.
  verify.logical.arg(normalize.rows)
  
  # SIMULATE DATA FROM MIXTURE
  # --------------------------
  # Argument "simtype" controls how the random numbers are generated.
  k <- floor(n/4)
  if (simtype == "n")
    x <- c(rnorm(n - 2*k),3*rnorm(k),6*rnorm(k))
  else if (simtype == "nt")
    x <- c(rnorm(n - 3*k),3*rnorm(k),5*rnorm(k),rt(k,df = 2))

  # SELECT VARIANCES FOR MIXTURE MODEL
  # ----------------------------------
  # Try to select a reasonable set of standard deviations that should
  # be used for the mixture model based on the values of x. This is
  # code is based on the autoselect.mixsd function from the ashr
  # package.
  smin <- 1/10
  if (all(x^2 < 1))
    smax <- 1
  else
    smax <- 2*sqrt(max(x^2 - 1))
  s <- c(0,logspace(smin,smax,m - 1))

  # CREATE LIKELIHOOD MATRIX
  # ------------------------
  # Entry (i,j) of the conditional likelihood matrix is equal to
  # N(0,se[i]^2 + s[j]^2), the normal density with zero mean and
  # variance se[i]^2 + s[j]^2, where se[i] is the standard error 
  # assigned to sample i. Here, all s.e.'s are assumed to be 1.
  L <- matrix(0,n,m)
  for (j in 1:m)
    L[,j] <- dnorm(x,sd = sqrt(1 + s[j]^2))

  # NORMALIZE LIKELIHOOD MATRIX
  # ---------------------------
  # Normalize the rows of the likelihood matrix so that the largest
  # entry in each row is 1.
  if (normalize.rows)
    L <- L / apply(L,1,max)

  # Return the simulated data points (x), the standard deviations
  # specifying the mixture prior (s), and the conditional likelihood
  # matrix (L).
  return(list(x = x,s = s,L = L))
}
