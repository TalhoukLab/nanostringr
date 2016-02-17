# Verifies that the input data have the correct column names
check_colnames <- function(...) {
  ind <- sapply(list(...), function(x) names(x)[1:3]) ==
    c("Code.Class", "Name", "Accession")
  if (all(ind))
    return(TRUE)
  else
    stop(paste(sapply(match.call()[-1], deparse)[
      which(!apply(ind, 2, function(x) all(x == TRUE)))],
      "does not have the correct column names.\n"))
}
