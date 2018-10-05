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
  rawdat <- raw.data[, -(1:3), drop = FALSE]
  rownames(rawdat) <- raw.data$Name
  hks <- raw.data$Code.Class == "Housekeeping"
  refs <- raw.data$Code.Class != "Endogenous"
  if (is.logged == FALSE) {
    rawdat <- rawdat + corr
    logHK <- apply(log2(rawdat[hks, , drop = FALSE]), 2, mean)
    logXpr <- log2(rawdat[!refs, , drop = FALSE])
  } else {
    logHK <- apply(rawdat[hks, , drop = FALSE], 2, mean)
    logXpr <- rawdat[!refs, , drop = FALSE]
  }
  if (ncol(logXpr) == 1) {
    norm <- logXpr - logHK
  } else {
    norm <- t(apply(logXpr, 1, function(x) x - logHK))
  }
  normdat <- cbind(raw.data[!refs, 1:3], norm)
  rownames(normdat) <- raw.data$Name[!refs]
  normdat
}
