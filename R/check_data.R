# Verifies that the input data for ratioMethod are valid
check_data <- function(...) {
  ind <- sapply(list(...), function(x) class(x) %in% c("matrix", "data.frame"))
  if (all(ind))
    return(TRUE)
  else
    stop(paste(sapply(match.call()[-1], deparse)[which(!ind)],
               "is not of class 'matrix' or 'data.frame'.\n"))
}
