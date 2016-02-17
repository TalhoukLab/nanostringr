# Verifies that all input data have the same number of columns
check_ncol <- function(...) {
  if (length(unique(sapply(list(...), ncol))) == 1)
    return(TRUE)
  else
    stop('All input data must have the same number of columns.')
}
