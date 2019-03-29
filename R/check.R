# Verifies that the input data have the correct column names
check_colnames <- function(...) {
  ind <- vapply(list(...), function(x) names(x)[1:3] ==
                  c("Code.Class", "Name", "Accession"), logical(3))
  if (all(ind)) {
    return(TRUE)
  } else{
    args <- vapply(match.call()[-1], deparse, character(1))
    stop(paste("\n", args[apply(ind, 2, function(x) !all(x))],
               "does not have the correct column names."))
  }
}

# Verifies that the input data for ratioMethod are valid
check_data <- function(...) {
  ind <- vapply(list(...), inherits, what = c("matrix", "data.frame"),
                logical(1))
  if (all(ind)) {
    return(TRUE)
  } else {
    args <- vapply(match.call()[-1], deparse, character(1))
    stop(paste("\n", args[!ind], "is not of class 'matrix' or 'data.frame'."))
  }
}

# Verifies that the Code.Class includes Endogenous and Housekeeping genes
check_genes <- function(x) {
  if (all(c("Endogenous", "Housekeeping") %in% x$Code.Class)) {
    return(TRUE)
  } else {
    stop("There are no Housekeeping genes in Code.Class.")
  }
}

# Verifies that all input data have the same number of columns
check_ncol <- function(...) {
  if (length(unique(vapply(list(...), ncol, integer(1)))) == 1) {
    return(TRUE)
  } else {
    stop('All input data must have the same number of columns.')
  }
}
