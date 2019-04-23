context("Parse RCC files")

rcc_file <- system.file("extdata", "example.RCC", package = "nanostringr")

test_that("reading from directory outputs list of counts and attributes", {
  rcc_dir <- dirname(rcc_file)
  rcc_list <- read_rcc(rcc_dir)
  expect_type(rcc_list, "list")
  expect_named(rcc_list, c("raw", "exp"))
  expect_length(rcc_list, 2)
})

test_that("counts data has right dimension", {
  counts <- parse_counts(rcc_file)
  expect_equal(nrow(counts), 354)
  expect_equal(ncol(counts), 4)
  expect_identical(
    head(names(counts), 3),
    c("Code.Class", "Name", "Accession")
  )
})

test_that("attributes data has right length and class", {
  attrs <- parse_attributes(rcc_file)
  expect_is(unlist(attrs[1:5]), "character")
  expect_is(unlist(attrs[6:8]), "numeric")
})
