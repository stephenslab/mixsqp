context("mixSQP")

test_that("Version number in mixSQP with verbose = TRUE is correct",{
  data(tacks)
  out <- capture.output(mixSQP(tacks$L,tacks$w))
  x   <- unlist(strsplit(out[1]," "))[[3]]
  expect_equal(packageDescription("mixSQP")$Version,x)
})

test_that(paste("mixSQP gives correct solutions for 2 x 2 and",
                "2 x 3 likelihood matrices"),{
  e <- 1e-8

  # In this first example, the correct solution is (1/2,1/2).
  L   <- rbind(c(1,e),
               c(e,1))
  capture.output(out <- mixSQP(L))
  expect_equal(out$x,c(0.5,0.5),tolerance = 1e-8)
  
  # In this second example, any solution of the form (x1,x2,0) gives
  # the same value for the objective.
  L    <- rbind(c(1,1,e),
                c(1,1,1))
  capture.output(out1 <- mixSQP(L,x0 = c(1,1,0)))
  capture.output(out2 <- mixSQP(L,x0 = c(0,1,1)))
  expect_equal(out1$status,mixsqp.status.converged)
  expect_equal(out2$status,mixsqp.status.converged)
  expect_equal(out1$x[3],0,tolerance = 1e-8)
  expect_equal(out2$x[3],0,tolerance = 1e-8)
})

test_that(paste("mixSQP and KWDual return the same solution for",
                "1000 x 10 likelihood matrix"),{

  # The REBayes package is required to run this test.
  skip_if_not_installed("REBayes")
                    
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
  capture.output(out2 <- mixSQP(L))

  # The outputted solutions, and the objective values at those
  # solutions, should be nearly identical. Also check that the
  # solution entries are labeled correctly.
  expect_equal(out2$status,mixsqp.status.converged)
  expect_equal(names(out1$x),colnames(L))
  expect_equal(names(out2$x),colnames(L))
  expect_equal(out1$x,out2$x,tolerance = 1e-6)
  expect_equal(out1$value,out2$value,tolerance = 1e-8)
})

test_that(paste("mixSQP & KWDual return the same solution for",
                "1000 x 10 likelihood matrix with unequal row weights"),{

  # The REBayes package is required to run this test.
  skip_if_not_installed("REBayes")
  
  # Simulate a 1000 x 10 likelihood matrix, and different weights for
  # the rows of this matrix.
  set.seed(1)
  L <- simulatemixdata(1000,10)$L
  w <- runif(1000)
  w <- w/sum(w)
  
  # Apply KWDual and mixSQP to the data set.
  out1 <- mixKWDual(L,w)
  capture.output(out2 <- mixSQP(L,w))

  # The outputted solutions, and the objective values at those
  # solutions, should be nearly identical.
  expect_equal(out2$status,mixsqp.status.converged)
  expect_equal(out1$x,out2$x,tolerance = 2e-8)
  expect_equal(out1$value,out2$value,tolerance = 2e-8)
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
  capture.output(out1 <- mixSQP(L1,w))
  capture.output(out2 <- mixSQP(L2,w))

  # The outputted solutions should be nearly identical (although the
  # values of the objectives will be different).
  expect_equal(out1$status,mixsqp.status.converged)
  expect_equal(out2$status,mixsqp.status.converged)
  expect_equal(out1$x,out2$x,tolerance = 1e-8)
})
                    
test_that(paste("mixSQP gives correct solution for Beckett & Diaconis",
                "tack rolling example"),{
  data(tacks)
  L <- tacks$L
  w <- tacks$w
  capture.output(out <- mixSQP(L,w))

  # The mixSQP solution should be very close to the REBayes solution
  # and, more importantly, the quality of the mixSQP solution should
  # be higher.
  expect_equal(out$status,mixsqp.status.converged)
  expect_equal(tacks$x,out$x,tolerance = 5e-4)
  expect_lte(out$value,mixobjective(L,tacks$x,w))
})

# This is mainly to test post-processing of the output when the
# algorithm reaches the maximum number of iterations. This example is
# used in one of the other tests above.
test_that("mixSQP does not report an error with convergence failure",{
  e <- 1e-8
  L <- rbind(c(1,1,e),
             c(1,1,1))
  capture_output(out <- mixSQP(L,x0 = c(0,1,1),maxiter.sqp = 3))
  expect_equal(out$status,mixsqp.status.didnotconverge)
  expect_equal(dim(out$data),c(3,6))
})

# This test comes from Issue #3.
test_that(paste("mixSQP gives correct solution for \"short and fat\" matrix,",
                "even when linear systems in active-set method are not",
                "necessarily s.p.d."),{
  set.seed(1)
  L    <- matrix(rgamma(1000,1,1),nrow = 10)
  out1 <- mixKWDual(L)
  capture.output(out2 <- mixSQP(L))
  capture.output(out3 <- mixSQP(L,delta = 0))
  
  # The mixSQP solution should be very close to the REBayes solution
  # and, more importantly, the quality of the mixSQP solution should
  # be very similar, even when the Newton search direction in the
  # active-set method is not necessarily unique (i.e., the Hessian is
  # not s.p.d.).
  expect_equal(out2$status,mixsqp.status.converged)
  expect_equal(out3$status,mixsqp.status.converged)
  expect_equal(out1$x,out2$x,tolerance = 1e-8)
  expect_equal(out1$x,out3$x,tolerance = 1e-8)
  expect_equal(out1$value,out2$value,tolerance = 1e-8)
  expect_equal(out1$value,out3$value,tolerance = 1e-8)
})

# This test comes from Issue #5.
test_that(paste("mixSQP converges, and outputs correct solution, for example",
                "in which the \"dual residual\" never reaches exactly zero"),{

  # Generate the data set for testing.
  set.seed(1)
  n    <- 1e5
  m    <- 12
  L    <- simulatemixdata(n,m)$L

  # Here we also check convergence for the case when no numerical
  # stability measure is used for the active-set linear systems (i.e.,
  # delta = 0).
  out1 <- mixKWDual(L)
  capture.output(out2 <- mixSQP(L,convtol.sqp = 0,maxiter.sqp = 20))
  capture.output(out3 <- mixSQP(L))
  capture.output(out4 <- mixSQP(L,delta = 0))
  
  # When the mixSQP iterates converge, they should be very close to
  # the REBayes solution, and when convtol.sqp = 0, the mixSQP
  # algorithm should report that it failed to converge in this
  # example.
  expect_equal(out2$status,"exceeded maximum number of iterations")
  expect_equal(out3$status,mixsqp.status.converged)
  expect_equal(out4$status,mixsqp.status.converged)
  expect_equal(out1$x,out3$x,tolerance = 1e-7)
  expect_equal(out1$x,out4$x,tolerance = 1e-7)
  expect_equal(out1$value,out3$value,tolerance = 1e-8)
})
