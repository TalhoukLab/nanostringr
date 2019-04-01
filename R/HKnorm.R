#' Normalization to Housekeeping Genes
#'
#' Normalizes the gene expression of NanoString nCounter data to housekeeping
#' genes. This is done by subtracting the average log housekeeping gene
#' expression from the expression level of every gene in each sample.
#'
#' @param raw.data matrix of raw counts obtained from nCounter (rows are genes).
#'   The first three columns must be labeled: c("Code.Class", "Name",
#'   "Accession") and contain that information.
#' @param is.logged logical; If `TRUE`, normalization has already been done
#'   on log base 2 scale, no need log the data
#' @param corr small correction to avoid error
#' @return A matrix of log normalized data in the same format but without
#'   reference genes.
#' @author Aline Talhouk
#' @export
#' @examples
#' HKnorm(ovd.r)
#' HKnorm(ovd.r, is.logged = TRUE)
HKnorm <- function(raw.data, is.logged = FALSE, corr = 0.0001) {
  assertthat::assert_that(check_colnames(raw.data))
  assertthat::assert_that(check_genes(raw.data))
  gxdat <- raw.data[raw.data$Code.Class == "Endogenous", -1:-3, drop = FALSE]
  hkdat <- raw.data[raw.data$Code.Class == "Housekeeping", -1:-3, drop = FALSE]
  if (is.logged) {
    hkmeans <- colMeans(hkdat)
  } else {
    gxdat <- log2(gxdat + corr)
    hkmeans <- colMeans(log2(hkdat + corr))
  }
  norm <- purrr::map2_df(gxdat, hkmeans, `-`)
  normdat <- raw.data[raw.data$Code.Class == "Endogenous", 1:3] %>%
    dplyr::bind_cols(norm)
  normdat
}
