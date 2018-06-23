#' @title testdata
#' 
#' @description samples L by means of nonparametric empirical Bayes approach for estimation of mixture proportions.
#' 
#' @param n nrow
#' @param m ncol
#' @param mix_type a type of mixture: "mix_n" for normal mixture, "mix_nt" for normal and t mixture.
#' 
#' @examples
#' n  <- 1e5
#' m  <- 20
#' L  <- testdata(n,m, mix_type = "mix_n") # Create some simulated data
#' 
#' @export
testdata <- function(n,m, mix_type = "mix_n"){
  
  temp = floor(n/4)
  if (mix_type == "mix_nt"){
    x = c(rnorm(n-3*temp),3*rnorm(temp),5*rnorm(temp),rt(temp,df = 2))
  } else if (mix_type == "mix_n"){
    x =  c(rnorm(n-2*temp),3*rnorm(temp),6*rnorm(temp))
  } else{
    stop("mixture type is not correctly defined")
  }
  
  smin = 1/10;
  if (all(x^2 < 1)){
    smax = 1
  } else{
    smax = min(sqrt(max(x^2 - 1)), 500); 
  }
  
  
  grid = c(0,10^seq(log10(smin),log10(smax),length = m-1))
  L <- matrix(0,n,m)
  
  for (j in 1:m)
    L[,j] <- dnorm(x,sd = sqrt(1 + grid[j]^2))
  
  L = L/apply(L,1,max)
  return(L)
}