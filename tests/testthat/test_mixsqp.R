context("mixSQP")

test_that("mixSQP and KWDual return the same solution",{

  # Apply KWDual and mixSQP to a 1,000 x 10 matrix.
  set.seed(1)
  L        <- simulatemixdata(1000,10)$L
  x.kwdual <- mixKWDual(L)$x
  x.mixsqp <- mixSQP(L,verbose = FALSE)$x
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
