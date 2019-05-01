#' Normalization to Housekeeping Genes
#'
#' Normalizes the gene expression of NanoString nCounter data to housekeeping
#' genes. This is done by subtracting the average log housekeeping gene
#' expression from the expression level of every gene in each sample.
#'
#' @inheritParams NanoStringQC
#' @param is.logged logical; If `TRUE`, normalization has already been done on
#'   log base 2 scale, no need log the data
#' @param corr small correction to avoid error
#' @return data frame of log normalized data in the same format but without
#'   reference genes
#' @author Aline Talhouk, Derek Chiu
#' @export
#' @examples
#' HKnorm(ovd.r)
#' HKnorm(ovd.r, is.logged = TRUE)
HKnorm <- function(raw, is.logged = FALSE, corr = 1e-04) {
  assertthat::assert_that(check_colnames(raw))
  assertthat::assert_that(check_genes(raw))
  gxdat <- raw[raw$Code.Class == "Endogenous", -1:-3, drop = FALSE]
  hkdat <- raw[raw$Code.Class == "Housekeeping", -1:-3, drop = FALSE]
  if (is.logged) {
    hkmeans <- colMeans(hkdat)
  } else {
    gxdat <- log2(gxdat + corr)
    hkmeans <- colMeans(log2(hkdat + corr))
  }
  norm <- purrr::map2_df(gxdat, hkmeans, `-`)
  normdat <- raw[raw$Code.Class == "Endogenous", 1:3] %>%
    dplyr::bind_cols(norm)
  normdat
}
