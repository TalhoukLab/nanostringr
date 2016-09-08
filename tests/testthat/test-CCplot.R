context("CCplot")

set.seed(12)
a1 <- rnorm(20)
a2 <- rnorm(20) + 2
b1 <- rnorm(19)

test_that("Measurements must be same length", {
  expect_error(CCplot(a1, a2), NA)
  expect_error(CCplot(a1, b1))
})

test_that("Outputs scatterplot or MAplot", {
  expect_error(CCplot(a1, a2, Ptype = "scatter"), NA)
  expect_error(CCplot(a1, a2, Ptype = "MAplot"), NA)
})

test_that("Metrics have correct structure", {
  m <- CCplot(a1, a2, metrics = TRUE)
  expect_is(m, "numeric")
  expect_length(m, 3)
})
