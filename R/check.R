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

# Verifies that the input data for ratioMethod are valid
check_data <- function(...) {
  ind <- sapply(list(...), function(x) class(x) %in% c("matrix", "data.frame"))
  if (all(ind))
    return(TRUE)
  else
    stop(paste(sapply(match.call()[-1], deparse)[which(!ind)],
               "is not of class 'matrix' or 'data.frame'.\n"))
}

# Verifies that the Code.Class includes Endogenous and Housekeeping genes
check_genes <- function(x) {
  if (all(c("Endogenous", "Housekeeping") %in% x$Code.Class))
    return(TRUE)
  else
    stop("There are no Housekeeping genes in Code.Class.")
}

# Verifies that all input data have the same number of columns
check_ncol <- function(...) {
  if (length(unique(sapply(list(...), ncol))) == 1)
    return(TRUE)
  else
    stop('All input data must have the same number of columns.')
}
