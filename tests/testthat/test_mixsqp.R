context("mixSQP")

test_that(paste("mixSQP gives correct solutions for 2 x 2 and",
                "2 x 3 likelihood matrices"),{
  e <- 1e-8

  # In this first example, the correct solution is (1/2,1/2).
  L   <- rbind(c(1,e),
               c(e,1))
  out <- mixSQP(L,verbose = FALSE)
  expect_equal(out$x,c(0.5,0.5),tolerance = 1e-8)
  
  # In this second example, any solution of the form (x1,x2,0) gives
  # the same value for the objective.
  L    <- rbind(c(1,1,e),
                c(1,1,1))
  out1 <- mixSQP(L,x0 = c(1,1,0),verbose = FALSE)
  out2 <- mixSQP(L,x0 = c(0,1,1),verbose = FALSE)
  expect_equal(out1$x[3],0,tolerance = 1e-8)
  expect_equal(out2$x[3],0,tolerance = 1e-8)
})

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

test_that(paste("mixSQP returns the same solution regardless whether",
                "the likelihood matrix is normalized"),{
  
  # Simulate two 100 x 10 likelihood matrices---one normalized and one
  # unnormalized---and different weights for the rows of this matrix.
  set.seed(1)
  L1 <- simulatemixdata(100,10,normalize.rows = TRUE)$L
  set.seed(1)
  L2 <- simulatemixdata(100,10,normalize.rows = FALSE)$L
  w  <- runif(100)
  w  <- w/sum(w)

  # Apply mixSQP to normalized and unnormalized data sets.
  out1 <- mixSQP(L1,w,verbose = FALSE)
  out2 <- mixSQP(L2,w,verbose = FALSE)

  # The outputted solutions should be nearly identical (although the
  # values of the objectives will be different).
  expect_equal(out1$x,out2$x,tolerance = 1e-8)
})
                    
