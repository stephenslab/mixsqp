# This file defines functions for simulating data sets.

# Simulate n random numbers X generated as follows: (1) with
# probability w[i], sample U from the univariate normal with zero mean
# and standard deviation s[i]; (2) sample X from the univariate normal
# with mean U and standard deviation se. This is the "Extreme
# Deconvolution" model for the special case when d = 1.
# 
# Input arguments w and and se must be vectors of the same length, and
# se specifies the standard devation in the "noise" for each of the
# samples, so it should be a vector of length n.
datasim.norm <- function (w, s, se) {

  # Get the number of mixture components (k) and the number of samples
  # to simulate (n).
  k <- length(w)
  n <- length(se)
    
  # Draw the source mixture component for each sample.
  i <- sample(k,n,replace = TRUE,prob = w)
  
  # Draw U for each sample.
  u <- rnorm(n = n,sd = s[i])

  # Draw X for each sample.
  return(rnorm(n = n,mean = u,sd = se))
}
