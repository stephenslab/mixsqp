context("mixSQP")

test_that("mixSQP and KWDual return the same solution",{

  # Apply KWDual and mixSQP to a 1,000 x 10 matrix.
  set.seed(1)
  L        <- simulatemixdata(1000,10)$L
  x.kwdual <- mixKWDual(L)$x
  x.mixsqp <- mixSQP(L,verbose = FALSE)$x
  expect_equal(x.kwdual,x.mixsqp,tolerance = 1e-5)
})
