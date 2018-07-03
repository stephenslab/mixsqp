context("mixSQP")

test_that("mixSQP and KWDual return the same solution",{

  # Simulate a 1,000 x 10 likelihood matrix. Note that I add row and
  # column names to the matrix to check that the column names are
  # retained in the solution vector.
  set.seed(1)
  n <- 1000
  m <- 10
  L <- simulatemixdata(n,m)$L
  rownames(L) <- paste0("x",1:n)
  colnames(L) <- paste0("s",1:m)
  
  # Apply KWDual and mixSQP to the data set.
  x.kwdual <- mixKWDual(L)$x
  x.mixsqp <- mixSQP(L,verbose = FALSE)$x
  expect_equal(x.kwdual,colnames(L))
  expect_equal(x.mixsqp,colnames(L))
  expect_equal(x.kwdual,x.mixsqp,tolerance = 1e-5)
})

test_that(paste("mixSQP & KWDual return the same solution for data",
                "with unequal weights"),{

  # Apply KWDual and mixSQP to a 100 x 10 matrix.
  set.seed(1)
  L        <- simulatemixdata(100,10)$L
  w        <- runif(100)
  w        <- w/sum(w) 
  out      <- REBayes::KWDual(L,rep(1,10),w)
  x <- out$f

  # Make sure the solution is (primal) feasible.
  x[x < 0] <- 0
  x        <- x/sum(x)

  # x.kwdual <- mixKWDual(L)$x
  x.mixsqp <- mixSQP(L,w = w,verbose = TRUE)$x
  expect_equal(x.kwdual,x.mixsqp,tolerance = 1e-5)
})
