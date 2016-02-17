#' Normalization to Housekeeping Genes (Single-Patient sample)
#'
#' Normalizes the gene expression of NanoString nCounter data to housekeeping genes.
#' This is done by subtracting the average log housekeeping gene expression.
#'
#' @param raw.data matrix of raw counts obtained from nCounter (rows are genes).
#' The first three columns must be labeled: c("Code.Class", "Name", "Accession")
#' and contain that information.
#' @param corr small correction on normalization
#' @return A matrix of log normalized data in the same format but without reference genes.
#' @author Aline Talhouk
#' @export
#' @examples
#' library(NanoStringNorm)
#' data(NanoString)
#' NanoString.mRNA[NanoString.mRNA$Name %in%
#' c('Eef1a1','Gapdh','Hprt1','Ppia','Sdha'), 'Code.Class'] <- 'Housekeeping'
#'
#' HKnorm(NanoString.mRNA)
HKnorm <- function(raw.data, corr = 0.0001) {
  assertthat::assert_that(check_colnames(raw.data))
  assertthat::assert_that(check_genes(raw.data))
  rawdat <- raw.data[, -(1:3)] + corr
  rownames(rawdat) <- raw.data$Name
  hks <- raw.data$Code.Class == "Housekeeping"
  refs <- raw.data$Code.Class != "Endogenous"
  logHK <- apply(log2(rawdat[hks, ]), 2, mean)
  logXpr <- log2(rawdat[!refs, ])
  norm <- t(apply(logXpr, 1, function(x) x - logHK))
  normdat <- cbind(raw.data[!refs, 1:3], norm)
  rownames(normdat) <- raw.data$Name[!refs]
  return(normdat)
}
