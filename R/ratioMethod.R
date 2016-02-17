#' Ratio method for batch adjustment
#'
#' Batch adjustment by taking the ratio of a reference sample
#' @param Y2 data run in the second batch, samples are rows and genes are columns
#' @param R1 reference data run in the first batch
#' @param R2 reference data run in the second batch
#' @return The Y2 data adjusted to batch 1
#' @author Aline Talhouk
#' @export
#' @examples
#' set.seed(12)
#' A <- matrix(rnorm(120), ncol = 10)
#' B <- matrix(rnorm(80), ncol = 10)
#' C <- matrix(rnorm(50), ncol = 10)
#' ratioMethod(A, B, C)
ratioMethod <- function(Y2, R1, R2) {
  assertthat::assert_that(check_data(Y2, R1, R2))
  assertthat::assert_that(check_ncol(Y2, R1, R2))
  m <- apply(R1, 2, mean) - apply(R2, 2, mean)
  Y2new <- t(apply(Y2, 1, function(x) x + m))
  return(Y2new)
}
