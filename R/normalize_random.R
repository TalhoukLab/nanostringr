#' Normalize data using random reference samples
#'
#' Normalize nanostring gene expression using randomly chosen samples
#' for the reference-based approach for batch adjustment.
#'
#' The number of randomly chosen numbers can be selected, and optionally a
#' `strata` can be specified such that `n` reference samples are selected from
#' each level (like a stratified bootstrap). In relation to the reference
#' method, the random samples removed from `ref` form `R1`, the random samples
#' removed from `x` form `R2`, and the remaining samples from `x` form `Y`. See
#' `refMethod()` for details.
#'
#' In subsequent analyses, we refer to a method using `normalize_random(n)` as
#' the "Random n" method.
#'
#' @param x target data
#' @param ref reference data
#' @param group_df grouping data used to select random reference samples
#' @param n number of random reference samples to select for normalization
#' @param strata character string of grouping variable found in `group_df` for
#'   stratified random sampling. If `strata` has `k` levels, then a total of `n * k`
#'   random samples are selected. The default uses no grouping variable and thus
#'   only simple random sampling is performed.
#' @param seed random seed for reproducibility
#' @return normalized gene expression
#'
#' @author Derek Chiu
#' @export
normalize_random <- function(x, ref, group_df, n = 1, strata = NULL, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (!is.null(strata)) {
    rand_df <- group_df %>%
      dplyr::distinct(.data$ottaID, .data[[strata]]) %>%
      dplyr::group_by(.data[[strata]]) %>%
      dplyr::slice_sample(n = n) %>%
      dplyr::ungroup()
  } else {
    rand_df <- group_df %>%
      dplyr::distinct(.data$ottaID) %>%
      dplyr::slice_sample(n = n)
  }

  ref_rand <- join_avg(ref, rand_df, "ottaID", "keep")
  x_rand <- join_avg(x, rand_df, "ottaID", "keep")
  x_counts <- join_avg(x, rand_df, "ottaID", "discard")
  x_norm <- as.data.frame(refMethod(x_counts, ref_rand, x_rand))
  return(x_norm)
}
