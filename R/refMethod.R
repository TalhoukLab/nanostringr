#' Reference-based approach for batch adjustment
#'
#' Batch adjustment by considering a measure relative to a reference sample
#'
#' @param Y data run in first or second batch, samples are rows and genes are
#'   columns. If correcting one batch only R1 is needed and would correspond to
#'   reference run in the same batch as Y, if calibrating one batch to the other
#'   Y represents the data from batch 2 and R1 would be reference run in batch 1
#'   and R2 would be reference from batch 2
#' @param R1 reference data run in the first batch
#' @param R2 reference data run in the second batch
#' @return The Y data adjusted calibrated to batch 1 (if two batches are
#'   presented) or the data with reference sample expression removed if only one
#'   data is provided
#' @author Aline Talhouk
#' @export
#' @examples
#' set.seed(12)
#' A <- matrix(rnorm(120), ncol = 10)
#' B <- matrix(rnorm(80), ncol = 10)
#' C <- matrix(rnorm(50), ncol = 10)
#' refMethod(A, B, C)
refMethod <- function(Y, R1, R2) {
  assertthat::assert_that(check_data(Y, R1, R2), check_ncol(Y, R1, R2))
  m <- colMeans(R1) - colMeans(R2)
  t(apply(Y, 1, `+`, m))
}
