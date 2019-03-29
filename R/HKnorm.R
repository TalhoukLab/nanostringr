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
  hkdat <- raw.data[raw.data$Code.Class == "Housekeeping", -1:-3, drop = FALSE]
  gxdat <- raw.data[raw.data$Code.Class == "Endogenous", -1:-3, drop = FALSE]
  if (!is.logged) {
    logHK <- colMeans(log2(hkdat + corr))
    logXpr <- log2(gxdat + corr)
  } else {
    logHK <- colMeans(hkdat)
    logXpr <- gxdat
  }
  if (ncol(logXpr) == 1) {
    norm <- logXpr - logHK
  } else {
    norm <- t(apply(logXpr, 1, function(x) x - logHK))
  }
  normdat <- cbind(raw.data[raw.data$Code.Class == "Endogenous", 1:3], norm)
  rownames(normdat) <- normdat$Name
  normdat
}
