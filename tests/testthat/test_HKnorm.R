context("HKnorm")

data(NanoString)

test_that("Error if no endogenous and housekeeping genes", {
  expect_error(HKnorm(NanoString.mRNA))
})

NanoString.mRNA[NanoString.mRNA$Name %in%
c('Eef1a1','Gapdh','Hprt1','Ppia','Sdha'), 'Code.Class'] <- 'Housekeeping'

test_that("Always returns data frame", {
  expect_is(HKnorm(NanoString.mRNA), "data.frame")
})
