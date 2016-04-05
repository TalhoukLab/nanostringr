context("NanoStringQC")

data("expQC", "hld.r", "hlo.r", "ovc.r", "ovd.r", "ovo.r")

test_that("Error if raw and exp different number of columns", {
  expect_error(NanoStringQC(hlo.r, expQC[expQC$geneRLF == "hlo.r", ]))
})

test_that("Always returns data frame", {
  expect_is(NanoStringQC(hlo.r, expQC[expQC$geneRLF == "HL1", ]), "data.frame")
})

test_that("ncol(raw) == nrow(exp) + 3", {
  expect_equal(ncol(hlo.r), nrow(expQC[expQC$geneRLF == "HL1", ]) + 3)
})
