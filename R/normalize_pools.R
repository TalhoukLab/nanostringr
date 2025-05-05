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
#' @param same_codeset logical; if `TRUE`, the target pool samples are collected
#'   from the same codset as the reference pool samples. This means the
#'   pool numbers from the file names use an older format where two-digit pools
#'   still pertain to pools 1, 2, or 3.
#' @return normalized gene expression
#'
#' @author Derek Chiu
#' @export
normalize_pools <- function(x, x_pools, ref_pools, p = 3, weigh = TRUE, same_codeset = FALSE) {
  pools <- paste0("Pool", seq_len(p))
  if (!same_codeset) {
    pool_nms <- rlang::set_names(pools, paste0(pools, "(?![0-9])"))
  } else {
    pool_nms <- rlang::set_names(pools)
  }
  if (weigh) {
    w <- pool_nms |>
      purrr::map_dbl(~ mean(grepl(., names(ref_pools), ignore.case = TRUE)))
  } else {
    w <- rlang::rep_named(pool_nms, 1 / length(pool_nms))
  }

  x_pools_mgx <- w |>
    purrr::imap(~ {
      .x * rowMeans(dplyr::select(x_pools, dplyr::matches(.y, perl = TRUE)))
    }) |>
    dplyr::bind_cols() |>
    rowSums() |>
    rlang::set_names(x_pools[["Name"]]) |>
    tibble::enframe(name = "Name", value = "x_exp")

  ref_pools_mgx <- dplyr::tibble(Name = rownames(ref_pools),
                                 ref_exp = unname(rowMeans(ref_pools)))

  x_norm <- x |>
    tidyr::pivot_longer(
      cols = dplyr::where(is.numeric),
      names_to = "FileName",
      values_to = "exp"
    ) |>
    dplyr::inner_join(x_pools_mgx, by = "Name") |>
    dplyr::inner_join(ref_pools_mgx, by = "Name") |>
    dplyr::mutate(norm_exp = exp + .data$ref_exp - .data$x_exp,
                  .keep = "unused") |>
    tidyr::pivot_wider(id_cols = "FileName",
                       names_from = "Name",
                       values_from = "norm_exp")
  return(x_norm)
}
