#' @title Create likelihood matrix from simulated data set
#' 
#' @description Simulate a data set, then compute the conditional
#'   likelihood matrix under a univariate normal likelihood and a
#'   mixture-of-normals prior. This models a simple nonparametric
#'   Empirical Bayes method applied to simulated data.
#' 
#' @param n Positive integer specifying the number of samples to
#'   generate and, consequently, the number of rows of the likelihood
#'   matrix L.
#' 
#' @param m Positiver integer specifying the number of mixture components.
#'
#' @param simtype The type of data to simulate. If \code{simtype =
#'   "n"}, simulate \code{n} random numbers from a mixture of three
#'   univariate normals with mean zero and standard deviation 1, 3 and
#'   6. If \code{simtype = "nt"}, simulate from a mixture of three
#'   univariate normals (with zero mean and standard deviations 1, 3 and
#'   5), and a t-distribution with 2 degrees of freedom.}
#'
#' @return 
#' 
#' @examples
#' n  <- 1e5
#' m  <- 20
#' L  <- testdata(n,m, mix_type = "mix_n") # Create some simulated data
#'
#' @importFrom stats rnorm
#' @importFrom stats dnorm
#' @importFrom stats rt
#' 
#' @export
#' 
simulatemixdata <- function (n, m, simtype = c("n","nt"), se = 1,
                             normalize.rows = TRUE) {

  # (1) CHECK INPUTS
  # ----------------
  # Input arguments n and m must be positive integers.
  if (n <= 0 | m <= 0)
    stop("Arguments \"n\" and \"m\" must be positive integers")
  n <- round(n)
  m <- round(m)

  # Get the choice of data to simulate.
  simtype <- match.arg(simtype)

  # TO DO: Check and process se.

  # TO DO: Check and process normalize.rows.
  
  # (2) SIMULATE DATA FROM MIXTURE
  # ------------------------------
  # 
  k <- floor(n/4)
  if (simtype == "n")
    x = c(rnorm(n - 2*k),3*rnorm(k),6*rnorm(k))
  else if (simtype == "nt")
    x = c(rnorm(n - 3*k),3*rnorm(k),5*rnorm(k),rt(k,df = 2))

  # (3) SELECT VARIANCES FOR MIXTURE MODEL
  # --------------------------------------
  # TO DO: Add more details here.
  smin = 1/10;
  if (all(x^2 < 1))
    smax = 1
  else{
    smax = min(sqrt(max(x^2 - 1)), 500); 
  }
  
  
  grid = c(0,10^seq(log10(smin),log10(smax),length = m-1))
  L <- matrix(0,n,m)
  
  for (j in 1:m)
    L[,j] <- dnorm(x,sd = sqrt(1 + grid[j]^2))
  
  L = L/apply(L,1,max)
  return(L)
}
