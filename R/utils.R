# magrittr placeholder
globalVariables(".")
globalVariables("where")

#' Join grouping variable to CodeSet
#'
#' Join grouping variable to CodeSet, then keep or discard the common ids, and
#' finally average the duplicate samples. Used as a helper function for
#' `normalize_random()`.
#'
#' @param cs nanostring expression of CodeSet
#' @param y grouping variable (e.g. histotype)
#' @param id identifier to join `cs` and `y` by
#' @param type the type of join to perform. If "keep", then a semi join keeps
#'   cases in `cs` that also appear in `y`. If "discard", then an anti join
#'   discards cases in `cs` that appear in `y`.
#' @author Derek Chiu
#' @noRd
join_avg <- function(cs, y, id, type = c("keep", "discard")) {
  type <- match.arg(type)
  join_fun <- switch(type,
                     keep = dplyr::semi_join,
                     discard = dplyr::anti_join)
  cs %>%
    join_fun(y, by = id) %>%
    dplyr::group_by(!!rlang::sym(id)) %>%
    dplyr::summarize_if(is.double, mean) %>%
    dplyr::ungroup() %>%
    tibble::column_to_rownames(id)
}
