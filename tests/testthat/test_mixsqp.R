context("mixSQP")

test_that(paste("mixSQP and KWDual return the same solution for",
                "1000 x 10 likelihood matrix"),{

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
  out1 <- mixKWDual(L)
  out2 <- mixSQP(L,verbose = FALSE)

  # The outputted solutions, and the objective values at those
  # solutions, should be nearly identical. Also check that the
  # solution entries are labeled correctly.
  expect_equal(names(out1$x),colnames(L))
  expect_equal(names(out2$x),colnames(L))
  expect_equal(out1$x,out2$x,tolerance = 1e-6)
  expect_equal(out1$value,out2$value,tolerance = 1e-8)
})

test_that(paste("mixSQP & KWDual return the same solution for",
                "1000 x 10 likelihood matrix with unequal row weights"),{

  # Simulate a 100 x 10 likelihood matrix, and different weights for
  # the rows of this matrix.
  set.seed(1)
  L <- simulatemixdata(100,10)$L
  w <- runif(100)
  w <- w/sum(w)
  
  # Apply KWDual and mixSQP to the data set.
  out1 <- mixKWDual(L,w)
  out2 <- mixSQP(L,w,verbose = FALSE)

  # The outputted solutions, and the objective values at those
  # solutions, should be nearly identical.
  expect_equal(out1$x,out2$x,tolerance = 1e-8)
  expect_equal(out1$value,out2$value,tolerance = 1e-8)
})
