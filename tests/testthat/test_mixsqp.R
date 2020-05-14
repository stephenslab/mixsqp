context("mixsqp")

# The Rmosek package on CRAN will not work with REBayes. This function
# is used for some of the tests to check whether the correct Rmosek
# package (the one downloaded from mosek.com) is installed.
skip_if_mixkwdual_doesnt_work <- function() {
  skip_if_not_installed("REBayes")
  skip_if_not_installed("Rmosek")
  skip_if(!is.element("mosek_lptoprob",getNamespaceExports("Rmosek")))
}

test_that("Version number in mixsqp with verbose = TRUE is correct",{
  L <- rbind(c(1,1,0),
             c(1,1,1))
  out <- capture.output(mixsqp(L))
  x   <- unlist(strsplit(out[1]," "))[4]
  expect_equal(packageDescription("mixsqp")$Version,x)
})

test_that(paste("mix-SQP allows zero likelihoods, but reports an error",
                "when initial estimate does not satisfy L*x > 0"),{
  L <- rbind(c(1,1,0),
             c(1,1,1))
  capture.output(out <- mixsqp(L))
  expect_equal(out$status,mixsqp:::mixsqp.status.converged)
  expect_error(mixsqp(L,x0 = c(0,0,1)))
})

test_that(paste("mix-SQP converges to correct solution even when initial",
                "estimate is very poor"),{
  e <- 1e-8
  L <- rbind(c(1,1,e),
             c(1,2,1))
  capture.output(out1 <- mixsqp(L))
  capture.output(out2 <- mixsqp(L,x0 = c(0,0,1)))
  expect_equal(out1$value,out2$value,tol = 1e-8)

  # This second example is particularly challenging because two of the
  # columns of the likelihood matrix are identical.
  L <- rbind(c(1,1,e),
             c(1,1,1))
  capture.output(out1 <- mixsqp(L))
  capture.output(out2 <- mixsqp(L,x0 = c(0,0,1)))
  expect_equal(out1$value,out2$value,tol = 1e-8)
})

test_that(paste("mix-SQP gives correct solutions for 2 x 2 and",
                "2 x 3 likelihood matrices"),{

  # In this first example, the correct solution is (1/2,1/2).
  e <- 1e-8
  L <- rbind(c(1,e),
             c(e,1))
  capture.output(out <- mixsqp(L))
  expect_equal(out$x,c(0.5,0.5),tolerance = 1e-8,scale = 1)
  
  # In this second example, any solution of the form (x1,x2,0) gives
  # the same value for the objective, and the third mixture weight
  # should be exactly zero.
  L <- rbind(c(1,1,e),
             c(1,1,1))
  capture.output(out1 <- mixsqp(L,x0 = c(1,1,0),control = list(eps = 0)))
  capture.output(out2 <- mixsqp(L,x0 = c(0,1,1),control = list(eps = 0)))
  expect_equal(out1$status,mixsqp:::mixsqp.status.converged)
  expect_equal(out2$status,mixsqp:::mixsqp.status.converged)
  expect_equal(out1$x[3],0,tolerance = 0)
  expect_equal(out2$x[3],0,tolerance = 0)
})

test_that(paste("mix-SQP and KWDual return the same solution for",
                "1000 x 10 likelihood matrix, and mix-SQP correctly",
                "estimates the nonzeros"),{

  # Simulate a 1,000 x 10 likelihood matrix. Note that I add row and
  # column names to the matrix to check that the column names are
  # retained in the solution vector.
  set.seed(1)
  n <- 1000
  m <- 10
  L <- simulatemixdata(n,m)$L
  rownames(L) <- paste0("x",1:n)
  colnames(L) <- paste0("s",1:m)
  
  # Apply mix-SQP solver to the data set. Check that the solution
  # entries are labeled correctly.
  capture.output(out1 <- mixsqp(L))
  expect_equal(out1$status,mixsqp:::mixsqp.status.converged)
  expect_equal(names(out1$x),colnames(L))

  # Apply KWDual solver to the data set. 
  skip_if_mixkwdual_doesnt_work()
  out2 <- mixkwdual(L)
  expect_equal(names(out2$x),colnames(L))

  # The outputted solutions, and the objective values at those
  # solutions, should be nearly identical.
  expect_equal(out1$x,out2$x,tolerance = 1e-4,scale = 1)
  expect_equal(out1$value,out2$value,tolerance = 1e-8,scale = 1)

  # The very small mixture weights estimated by KWDual are all smaller
  # in the mix-SQP output.
  i <- which(out2$x < 0.001)
  expect_equivalent(out1$x[i] < out2$x[i],rep(TRUE,length(i)))
})

test_that(paste("mix-SQP & KWDual return the same solution for",
                "1000 x 10 likelihood matrix with unequal row weights"),{

  # Simulate a 1000 x 10 likelihood matrix, and different weights for
  # the rows of this matrix.
  set.seed(1)
  L <- simulatemixdata(1000,10)$L
  w <- runif(1000)
  w <- w/sum(w)
  
  # Apply mix-SQP solver to the data set.
  capture.output(out1 <- mixsqp(L,w))
  expect_equal(out1$status,mixsqp:::mixsqp.status.converged)

  # Apply KWDual solver to the data set.
  skip_if_mixkwdual_doesnt_work()
  out2 <- mixkwdual(L,w)
  
  # The outputted solutions, and the objective values at those
  # solutions, should be nearly identical.
  expect_equal(out1$x,out2$x,tolerance = 1e-4,scale = 1)
  expect_equal(out1$value,out2$value,tolerance = 1e-6,scale = 1)
})

test_that(paste("mix-SQP returns the same solution regardless whether",
                "the likelihood matrix is normalized"),{
  
  # Simulate two 100 x 10 likelihood matrices---one normalized and one
  # unnormalized---and different weights for the rows of this matrix.
  set.seed(1)
  L1 <- simulatemixdata(100,10,normalize.rows = TRUE)$L
  set.seed(1)
  L2 <- simulatemixdata(100,10,normalize.rows = FALSE)$L
  w  <- runif(100)
  w  <- w/sum(w)

  # Apply mix-SQP to normalized and unnormalized data sets.
  capture.output(out1 <- mixsqp(L1,w))
  capture.output(out2 <- mixsqp(L2,w))

  # The outputted solutions should be nearly identical (although the
  # values of the objectives will be different).
  expect_equal(out1$status,mixsqp:::mixsqp.status.converged)
  expect_equal(out2$status,mixsqp:::mixsqp.status.converged)
  expect_equal(out1$x,out2$x,tolerance = 1e-8,scale = 1)
})

test_that(paste("mix-SQP returns the correct solution when log-likelihoods",
                "are provided"),{
  
  # Simulate a 100 x 10 likelihood matrix as well as different weights
  # for the rows of this matrix.
  set.seed(1)
  L <- simulatemixdata(100,10,normalize.rows = FALSE)$L
  w <- runif(100)
  w <- w/sum(w)

  # Compute the log-likelihood matrix for the same data.
  set.seed(1)
  logL <- simulatemixdata(100,10,log = TRUE)$L
  
  # Apply mix-SQP to the likelihoods and the log-likelihoods.
  capture.output(out1 <- mixsqp(L,w))
  capture.output(out2 <- mixsqp(logL,w,log = TRUE))

  # The outputted solutions should be nearly identical (though the
  # values of the objectives will be different).
  expect_equal(out1$status,mixsqp:::mixsqp.status.converged)
  expect_equal(out2$status,mixsqp:::mixsqp.status.converged)
  expect_equal(out1$x,out2$x,tolerance = 1e-8,scale = 1)
})

test_that(paste("mix-SQP returns the same solution when using full data",
                "matrix and the low-rank SVD approximation"),{

  # Simulate a 1,000 x 80 likelihood matrix.
  set.seed(1)
  L <- simulatemixdata(1000,80)$L
  w <- runif(1000)
  w <- w/sum(w)

  # Apply mix-SQP to normalized and unnormalized data sets.
  capture.output(out1 <- mixsqp(L,w,control = list(tol.svd = 0)))
  capture.output(out2 <- mixsqp(L,w,control = list(tol.svd = 1e-8)))

  # The outputted solutions---and the gradients and Hessians at the
  # estimated solutions---should be nearly identical.
  expect_equal(out1$value,out2$value,tolerance = 1e-6,scale = 1)
  expect_equal(out1$x,out2$x,tolerance = 1e-6,scale = 1)
  expect_equal(out1$grad,out2$grad,tolerance = 1e-6,scale = 1)
  expect_equal(out1$hessian,out2$hessian,tolerance = 1e-6,scale = 1)
})

test_that("mix-SQP successfully \"escapes\" a sparse initial estimate",{
  set.seed(1)
  n <- 100
  m <- 10
  out <- simulate_data_koenker(n,m)
  L   <- out$L
  w   <- out$w
  x0  <- c(1,rep(0,m - 1))
  capture.output(fit1 <- mixsqp(L,w,x0,control = list(numiter.em = 0)))
  expect_equal(fit1$status,mixsqp:::mixsqp.status.converged)

  # Compare mix-SQP solution to KWDual.
  skip_if_mixkwdual_doesnt_work()
  fit2 <- mixkwdual(L,w)
  expect_equal(fit1$x,fit2$x,tolerance = 1e-4,scale = 1)
  expect_lte(fit1$value,fit2$value)
})

test_that(paste("mix-SQP gives correct solution for Beckett & Diaconis",
                "tack rolling example"),{

  data(tacks)
  L <- tacks$L
  w <- tacks$w
  capture.output(out <- mixsqp(L,w))

  # The objective value at the mix-SQP solution should be very close
  # to the objective value at the KWDual solution.
  expect_equal(out$status,mixsqp:::mixsqp.status.converged)
  expect_equal(mixobjective(L,out$x,w),mixobjective(L,tacks$x,w),
               tolerance = 1e-4)
})

# This is mainly to test post-processing of the output when the
# algorithm reaches the maximum number of iterations. This example is
# used in one of the other tests above.
test_that("mix-SQP gives warning with convergence failure",{
  e <- 1e-8
  L <- rbind(c(1,1,e),
             c(1,1,1))
  expect_warning(capture_output(out <- mixsqp(L,x0 = c(0,1,1),
                                  control = list(numiter.em  = 0,
                                                 maxiter.sqp = 1))))
  expect_equal(out$status,mixsqp.status.didnotconverge)
  expect_equal(dim(out$progress),c(1,7))
})

# This test comes from Issue #3.
test_that(paste("mix-SQP gives correct solution for \"short and fat\" matrix,",
                "even when linear systems in active-set method are not",
                "necessarily s.p.d."),{
  set.seed(1)
  L <- matrix(rgamma(1000,1,1),nrow = 10)
  capture.output(out1 <- mixsqp(L))

  # The mix-SQP solution should be very close to the KWDual solution
  # and, more importantly, the quality of the mix-SQP solution should
  # be very similar, even when the Newton search direction in the
  # active-set method is not necessarily unique (i.e., the Hessian is
  # not s.p.d.).
  skip_if_mixkwdual_doesnt_work()
  out2 <- mixkwdual(L)
  expect_equal(out1$x,out2$x,tolerance = 1e-4,scale = 1)
  expect_equal(out1$value,out2$value,tolerance = 1e-6,scale = 1)
})

# This test comes from Issue #5.
test_that(paste("mix-SQP converges, and outputs correct solution, for example",
                "in which the \"dual residual\" never reaches exactly zero"),{

  # Generate the data set for testing.
  set.seed(1)
  n <- 1e5
  m <- 12
  L <- simulatemixdata(n,m)$L

  # Here we also check convergence for the case when no numerical
  # stability measure is used for the active-set linear systems (i.e.,
  # eps = 0). Also, when convtol.sqp = 0, the mix-SQP algorithm should
  # report that it failed to converge in this example.
  suppressWarnings(capture.output(out1 <-
    mixsqp(L,control = list(eps = 0,convtol.sqp = 0,maxiter.sqp = 10))))
  capture.output(out2 <- mixsqp(L))
  capture.output(out3 <- mixsqp(L,control = list(eps = 0)))
  expect_equal(out1$status,"exceeded maximum number of iterations")
  expect_equal(out2$status,mixsqp:::mixsqp.status.converged)
  expect_equal(out3$status,mixsqp:::mixsqp.status.converged)
  
  # When the mix-SQP iterates converge, they should be very close to
  # the KWDual solution.
  skip_if_mixkwdual_doesnt_work()
  out4 <- mixkwdual(L)
  expect_equal(out2$x,out4$x,tolerance = 1e-5,scale = 1)
  expect_equal(out3$x,out4$x,tolerance = 1e-5,scale = 1)
  expect_equal(out2$value,out4$value,tolerance = 1e-6,scale = 1)
  expect_equal(out3$value,out4$value,tolerance = 1e-6,scale = 1)
})

test_that(paste("Case is properly handled in which all columns except",
                "one are filled with zeros"),{
  set.seed(1)
  n       <- 200
  m       <- 10
  i       <- 7
  L       <- simulatemixdata(n,m)$L
  L[,-i]  <- 0
  xsol    <- rep(0,m)
  xsol[i] <- 1
  expect_warning(capture.output(out1 <- mixsqp(L)))
  expect_equal(out1$status,mixsqp:::mixsqp.status.didnotrun)
  expect_null(out1$progress)
  expect_equal(out1$x,xsol)

  skip_if_mixkwdual_doesnt_work()
  expect_warning(capture.output(out2 <- mixkwdual(L)))
  expect_equal(out2$x,xsol)
})

test_that("Case is properly handled in which one column of L is all zeros",{
  set.seed(1)
  n <- 200
  m <- 10
  i <- 7
  L <- simulatemixdata(n,m)$L

  # Run the mix-SQP algorithm when all columns have nonzeros.
  L[,i] <- 1e-8
  capture.output(out1 <- mixsqp(L))

  # Set one of the columns to be all zeros, and re-run the mix-SQP
  # algorithm.
  L[,i] <- 0
  expect_warning(capture.output(out2 <- mixsqp(L)))

  # The two solutions should be pretty much the same.
  expect_equal(out1$x,out2$x,tolerance = 1e-4,scale = 1)
  expect_equal(out1$value,out2$value,tolerance = 1e-8,scale = 1)

  # Also check KWDual solution.
  skip_if_mixkwdual_doesnt_work()
  expect_warning(out3 <- mixkwdual(L))
  expect_equal(out1$x,out3$x,tolerance = 0.001,scale = 1)
  expect_equal(out1$value,out3$value,tolerance = 1e-6,scale = 1)
})

test_that("Case is properly handled when L has only one column",{
  set.seed(1)
  L    <- matrix(runif(100))
  suppressWarnings(out1 <- mixsqp(L))
  suppressWarnings(out2 <- mixkwdual(L))
  expect_equal(out1$x,1)
  expect_equal(out2$x,1)
  expect_equal(out1$value,out2$value)
})

test_that(paste("mix-SQP converges to solution for \"flat\" objective even",
                "if initial progress is poor"),{
  set.seed(1)
  n <- 10000
  m <- 20
  L <- matrix(runif(n*m),n,m)
  capture.output(out <- mixsqp(L))
  expect_equal(out$status,mixsqp:::mixsqp.status.converged)
})

# This test comes from Issue #19.
test_that("mix-SQP converges in a more difficult example",{
  load("flashr.example.RData")
  capture.output(out1 <- mixsqp(L))
  expect_equal(out1$status,mixsqp:::mixsqp.status.converged)
  skip_if_mixkwdual_doesnt_work()
  out2 <- mixkwdual(L)
  expect_equal(out1$x,out2$x,tolerance = 0.001,scale = 1)
})

# This test comes from Issue #29.
test_that("mix-SQP converges despite poor initialization",{
  load("ashr.example.RData")
  capture.output(out1 <- mixsqp(L,x0 = x0))
  expect_equal(out1$status,mixsqp:::mixsqp.status.converged)
})

# This test comes from Issue #30.
test_that("mix-SQP works for difficult smashr example",{
  load("smashr.example.RData")
  capture.output(out1 <- mixsqp(L,w,x0))
  expect_equal(out1$status,mixsqp:::mixsqp.status.converged)
  skip_if_mixkwdual_doesnt_work()
  out2 <- mixkwdual(L)
  expect_equal(out1$x,out2$x,tolerance = 1e-6,scale = 1)
})

# This example was generated by running:
#
#   library(ashr)
#   set.seed(1)
#   logitp <- c(rep(0,800), runif(200,-3,3))
#   p      <- 1/(1 + exp(-logitp))
#   n      <- rep(100,1000)
#   x      <- rbinom(1000,n,p)
#   out    <- ash(rep(0,length(x)),1,lik = lik_binom(x,n,link = "logit"),
#                 mode = "estimate")
#
test_that("mix-SQP works for difficult ashr example with binomial likelihood",{
  skip_if_mixkwdual_doesnt_work()
  load("ashr.binom.example.RData")
  capture.output(out1 <- mixsqp(L))
  out2 <- mixkwdual(L)
  expect_equal(out1$value,out2$value,tolerance = 1e-8,scale = 1)
  expect_equal(out1$x,out2$x,tolerance = 1e-4,scale = 1)
})

# This test comes from stephens999/ashr Issue #76.
test_that("mix-SQP works for a difficult mashr example",{
  skip_if_mixkwdual_doesnt_work()
  skip_if_not(file.exists("mashr.example.RData"))
  load("mashr.example.RData")
  capture.output(out1 <- mixsqp(L,w))
  out2 <- mixkwdual(L,w)
  expect_equal(out1$value,out2$value,tolerance = 1e-8,scale = 1)
  expect_equal(out1$x,out2$x,tolerance = 1e-4,scale = 1)
})

# This test comes from Issue #42.
test_that(paste("mix-SQP converges, and gives solution that is equally as",
                "good as KWDual for difficult NPMLE problem"),{
  skip_if_not(file.exists("npmle.RData"))
  load("npmle.RData")
  capture.output(out1 <- mixsqp(L))
  expect_equal(out1$status,mixsqp:::mixsqp.status.converged)
  skip_if_mixkwdual_doesnt_work()
  out2 <- mixkwdual(L)

  # The likelihood surface is very "flat", so the solutions are not
  # expected to be the same; however, the quality of the mix-SQP
  # solution should be as good, or nearly as good, as the solution
  # returned by the KWDual solver.
  expect_equal(out1$value,out2$value,tolerance = 1e-4,scale = 1)
})
