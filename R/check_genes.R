# Verifies that the Code.Class includes Endogenous and Housekeeping genes
check_genes <- function(x) {
  if (all(c("Endogenous", "Housekeeping") %in% x$Code.Class))
    return(TRUE)
  else
    stop("There are no Housekeeping genes in Code.Class.")
}
