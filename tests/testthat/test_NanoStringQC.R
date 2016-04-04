context("NanoStringQC")

library(otta)
data(rawOVCA2, rawPROT, rawOTTA, annot)
cs1 <- rawOVCA2
cs2 <- rawPROT
cs3 <- rawOTTA
exp0 <- annot

test_that("Error if raw and exp different number of columns", {
  expect_error(NanoStringQC(cs1, exp0[exp0$geneRLF == "CS1", ]))
  expect_error(NanoStringQC(cs2, exp0[exp0$geneRLF == "CS2", ]))
  expect_error(NanoStringQC(cs3, exp0[exp0$geneRLF == "CS3", ]))
})

exp0$geneRLF <- as.character(factor(
  exp0$geneRLF, labels = c("HL1", "HL2", "HL3", "HuRef",
                           "CS3", "mini", "CS1", "CS2")))

test_that("Always returns data frame", {
  expect_is(NanoStringQC(cs1, exp0[exp0$geneRLF == "CS1", ]), "data.frame")
  expect_is(NanoStringQC(cs2, exp0[exp0$geneRLF == "CS2", ]), "data.frame")
  expect_is(NanoStringQC(cs3, exp0[exp0$geneRLF == "CS3", ]), "data.frame")
})

test_that("ncol(raw) == nrow(exp) + 3", {
  expect_equal(ncol(cs1), nrow(exp0[exp0$geneRLF == "CS1", ]) + 3)
  expect_equal(ncol(cs2), nrow(exp0[exp0$geneRLF == "CS2", ]) + 3)
  expect_equal(ncol(cs3), nrow(exp0[exp0$geneRLF == "CS3", ]) + 3)
})
