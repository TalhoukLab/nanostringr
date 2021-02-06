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

test_that("No errors if column names start with a number", {
  colnames(hlo.r) <- gsub("^X", "",  colnames(hlo.r))
  expect_error(NanoStringQC(hlo.r, expQC[expQC$geneRLF == "HL1", ]), NA)
})

test_that("Error if positive controls don't have concentrations in brackets", {
  hlo.r$Name[hlo.r$Code.Class == "Positive"] <-
    gsub("\\([0-9].*", "\\1", hlo.r$Name[hlo.r$Code.Class == "Positive"])
  expect_error(NanoStringQC(hlo.r, expQC[expQC$geneRLF == "HL1", ]))
})

test_that("Error if non-standard column names", {
  hlo.r.2 <- hlo.r
  names(hlo.r.2)[1:2] <- c("A", "B")
  expect_error(NanoStringQC(hlo.r.2, expQC[expQC$geneRLF == "HL1", ]))
})

test_that("Error if missing any gene types", {
  hlo.r.3 <- dplyr::filter(hlo.r, Code.Class != "Housekeeping")
  expect_error(NanoStringQC(hlo.r.3, expQC[expQC$geneRLF == "HL1", ]))

  hlo.r.4 <- dplyr::filter(hlo.r, Code.Class != "Endogenous")
  expect_error(NanoStringQC(hlo.r.4, expQC[expQC$geneRLF == "HL1", ]))

  hlo.r.5 <- dplyr::filter(hlo.r, Code.Class != "Negative")
  expect_error(NanoStringQC(hlo.r.5, expQC[expQC$geneRLF == "HL1", ]))

  hlo.r.6 <- dplyr::filter(hlo.r, Code.Class != "Positive")
  expect_error(NanoStringQC(hlo.r.6, expQC[expQC$geneRLF == "HL1", ]))
})
