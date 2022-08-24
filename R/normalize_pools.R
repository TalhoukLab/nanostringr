#' Normalize data using common pools
#'
#' Normalize nanostring gene expression using common pools between two CodeSets.
#'
#' The target and reference expression samples, as well the target and reference
#' pool samples all need to be specified. We recommend reweighing the target
#' pool samples when calculating the average expression by the distribution of
#' reference pools.
#'
#' @inheritParams normalize_random
#' @param x_pools target pool samples
#' @param ref_pools reference pool samples
#' @param p number of pool sample sets. Defaults to 3.
#' @param weigh logical; if `TRUE`, the average expression in `x_pools` is
#'   reweighed by the distribution of the `p` pool sample sets in `ref_pools`.
#' @return normalized gene expression
#'
#' @author Derek Chiu
#' @export
normalize_pools <- function(x, ref, x_pools, ref_pools, p = 3, weigh = TRUE) {
  pool_nms <- rlang::set_names(paste0("Pool", seq_len(p)))
  if (weigh) {
    w <- pool_nms %>%
      purrr::map_dbl(~ mean(grepl(., names(ref_pools), ignore.case = TRUE)))
  } else {
    w <- pool_nms %>%
      rlang::rep_named(1 / length(.))
  }

  x_pools_mgx <- w %>%
    purrr::imap_dfc(~ .x * rowMeans(dplyr::select(x_pools, dplyr::matches(.y)))) %>%
    dplyr::transmute(Name = x_pools[["Name"]], x_exp = rowSums(.))
  ref_pools_mgx <- dplyr::tibble(Name = rownames(ref_pools),
                                 ref_exp = unname(rowMeans(ref_pools)))

  x_val <- x %>%
    dplyr::select(.data$Name, setdiff(names(.), names(x_pools)))
  x_norm <- x_val %>%
    tidyr::pivot_longer(cols = where(is.numeric),
                        names_to = "FileName",
                        values_to = "exp") %>%
    dplyr::inner_join(x_pools_mgx, by = "Name") %>%
    dplyr::inner_join(ref_pools_mgx, by = "Name") %>%
    dplyr::mutate(be = .data$x_exp - .data$ref_exp) %>%
    dplyr::transmute(
      Name = factor(.data$Name, levels = unique(.data$Name)),
      .data$FileName,
      exp = .data$be + exp
    ) %>%
    tidyr::pivot_wider(names_from = "Name", values_from = "exp")
  return(x_norm)
}
