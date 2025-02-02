set.seed(1)
logXpr <- replicate(5, rnorm(6))
dim(logXpr)

logHK <- apply(logXpr[, 4:5], 1, mean)
logref <- apply(logXpr[5:6, ], 2, mean)

normtoHK <- sweep(sweep(logXpr, 1, logHK), 2, logref)
normtoPool <- sweep(sweep(logXpr, 2, logref), 1, logHK)

test_that("Same Order", {
  expect_equal(normtoHK, normtoPool)
})
