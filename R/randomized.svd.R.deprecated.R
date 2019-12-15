randomized.svd <- function(A, k=NULL, nu=NULL, nv=NULL, p=0, q=0, sdist="normal") {
  
  #Dim of input matrix
  m <- nrow(A)
  n <- ncol(A)
  
  #Flipp matrix, if wide
  if(m < n){
    A <- H(A)
    m <- nrow(A)
    n <- ncol(A)
    flipped <- TRUE
  } else flipped <- FALSE
  
  #Set target rank
  if(is.null(k)) k = n
  if(k > n) k <- n
  if(is.character(k)) stop("Target rank is not valid!")
  if(k < 1) stop("Target rank is not valid!")
  
  #Set oversampling parameter
  l <- round(k) + round(p)
  if(l > n) l <- n
  if(l < 1) stop("Target rank is not valid!")
  
  #Check if array is real or complex
  if(is.complex(A)) {
    isreal <- FALSE
  } else {
    isreal <- TRUE
  }
  
  #Set number of singular vectors
  if(is.null(nu)) nu <- k
  if(is.null(nv)) nv <- k
  if(nu < 0) nu <- 0
  if(nv < 0) nv <- 0
  if(nu > k) nu <- k
  if(nv > k) nv <- k
  if(flipped==TRUE) {
    temp <- nu
    nu <- nv
    nv <- temp
  }
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Generate a random sampling matrix O
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  O <- switch(sdist,
              normal = matrix(stats::rnorm(l*n), n, l),
              unif = matrix(stats::runif(l*n), n, l),
              rademacher = matrix(sample(c(-1,1), (l*n), replace = TRUE, prob = c(0.5,0.5)), n, l),
              stop("Selected sampling distribution is not supported!"))
  
  if(isreal==FALSE) {
    O <- O + switch(sdist,
                    normal = 1i * matrix(stats::rnorm(l*n), n, l),
                    unif = 1i * matrix(stats::runif(l*n), n, l),
                    rademacher = 1i * matrix(sample(c(-1,1), (l*n), replace = TRUE, prob = c(0.5,0.5)), n, l),
                    stop("Selected sampling distribution is not supported!"))
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Build sample matrix Y : Y = A * O
  #Note: Y should approximate the range of A
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Y <- A %*% O
  remove(O)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Orthogonalize Y using economic QR decomposition: Y=QR
  #If q > 0 perfrom q subspace iterations
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if( q > 0 ) {
    for( i in 1:q) {
      Y <- qr.Q( qr(Y, complete = FALSE) , complete = FALSE )
      Z <- crossprod(A , Y)
      Z <- qr.Q( qr(Z, complete = FALSE) , complete = FALSE )
      Y <- A %*% Z
    }#End for
    remove(Z)
  }#End if
  
  Q <- qr.Q( qr(Y, complete = FALSE) , complete = FALSE )
  remove(Y)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Project the data matrix a into a lower dimensional subspace
  #B := Q.T * A
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  B <- crossprod(Q , A)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Singular Value Decomposition
  #Note: B =: U * S * Vt
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  rsvdObj <- svd(B, nu=nu, nv=nv) # Compute SVD
  rsvdObj$d <- rsvdObj$d[1:k] # Truncate singular values
  
  if(nu != 0) rsvdObj$u <- Q %*% rsvdObj$u # Recover left singular vectors
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Flip SVD back
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(nu == 0){ rsvdObj$u <- NULL}
  if(nv == 0){ rsvdObj$v <- NULL}
  if(flipped == TRUE) {
    u_temp <- rsvdObj$u
    rsvdObj$u <- rsvdObj$v
    rsvdObj$v <- u_temp
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Return
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  class(rsvdObj) <- "rsvd"
  return(rsvdObj) 
  
} # End rsvd