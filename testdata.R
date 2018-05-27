testdata <- function(n,m, type = 1){
  
  temp = floor(n/4)
  if (type == 1){
    x = c(rnorm(n-3*temp),3*rnorm(temp),5*rnorm(temp),rt(temp,df = 2))
  } else{
    x =  c(rnorm(n-3*temp),3*rnorm(temp),6*rnorm(2*temp))
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