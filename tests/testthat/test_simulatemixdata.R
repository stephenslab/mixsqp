context("mixSQP")

test_that("simulatemixdata works",{
  dat <- simulatemixdata(1000,20)
  expect_equal(dim(dat$L),c(1000,20))
})
