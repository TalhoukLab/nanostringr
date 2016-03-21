#' Concordance Correlation Plot
#'
#' Plotting function for reliability measure.
#' @param method1 measurements obtained in batch 1 or using method 1
#' @param method2 measurements obtained in batch 2 or using method 2
#' @param Ptype type of plot to be outputted c("scatter", "MAplot")
#' @param metrics if return metrics is set to \code{TRUE}, returns Rc, Ca and R2
#' @param xlabel label for x axis
#' @param ylabel label for y axis
#' @param title title for the main plot
#' @param subtitle subtitle of plot
#' @param xrange range of x axis
#' @param yrange range of y axis
#' @param MArange MA range
#' @return Either a scatterplot or MA plot showing concordance correlation.
#' @author Aline Talhouk
#' @export
#' @examples
#' # Simulate normally distributed data
#' set.seed(12)
#' a1 <- rnorm(20) + 2
#' a2 <- a1 + rnorm(20, 0, 0.15)
#' a3 <- a1 + rnorm(20, 0, 0.15) + 1.4
#' a4 <- 1.5 * a1 + rnorm(20, 0, 0.15)
#' a5 <- 1.3 * a1 + rnorm(20, 0, 0.15) + 1
#' a6 <- a1 + rnorm(20, 0, 0.8)
#' par(mfrow = c(3, 2), mar = c(5.1, 4.1, 1.5, 1.5))
#'
#' # Scatterplots
#' CCplot(a1, a1, Ptype = "scatter", "X", "Y", "Perfect Agreement", subtitle = letters[1])
#' CCplot(a1, a2, Ptype = "scatter", "X", "Y", "Very Good Agreement", subtitle = letters[2])
#' CCplot(a1, a3, Ptype = "scatter", "X", "Y", "Location Shift", subtitle = letters[3])
#' CCplot(a1, a4, Ptype = "scatter", "X", "Y", "Scale Shift", subtitle = letters[4])
#' CCplot(a1, a5, Ptype = "scatter", "X", "Y", "Location and Scale Shift", subtitle = letters[5])
#' CCplot(a1, a6, Ptype = "scatter", "X", "Y", "Measurement Error", subtitle = letters[6])
#'
#' # MAplots
#' CCplot(a1, a1, Ptype = "MAplot", "X", "Y", "Perfect Agreement", subtitle = letters[1])
#' CCplot(a1, a2, Ptype = "MAplot", "X", "Y", "Very Good Agreement", subtitle = letters[2])
#' CCplot(a1, a3, Ptype = "MAplot", "X", "Y", "Location Shift", subtitle = letters[3])
#' CCplot(a1, a4, Ptype = "MAplot", "X", "Y", "Scale Shift", subtitle = letters[4])
#' CCplot(a1, a5, Ptype = "MAplot", "X", "Y", "Location and Scale Shift", subtitle = letters[5])
#' CCplot(a1, a6, Ptype = "MAplot", "X", "Y", "Measurement Error", subtitle = letters[6])
CCplot <- function(method1, method2, Ptype = "None", metrics = FALSE,
                   xlabel = "", ylabel = "", title = "", subtitle = "",
                   xrange = NULL, yrange = NULL, MArange = c(-3.5, 5.5)) {
  assertthat::assert_that(length(method1) == length(method2))
  tmp.ccc <- epiR::epi.ccc(method1, method2, ci = "z-transform",
                           conf.level = 0.95)
  cclab <- paste0("Rc: ", round(tmp.ccc$rho.c[,1], digits = 2), " (",
                 round(tmp.ccc$rho.c[, 2], digits = 2), " - ",
                 round(tmp.ccc$rho.c[, 3], digits = 2), ")")
  r2lab <- bquote(R : .( round(cor(method1, method2),2)))
  r2rob <- bquote(rM : .( round(ccaPP::corM(method1, method2),2)))
  Acc <- paste("Ca:", round(tmp.ccc$C.b, 2))
  loc <- paste0("Location shift:", round(tmp.ccc$l.shift, 2))
  scl <- paste0("Scale shift:", round(tmp.ccc$s.shift, 2))
  z <- lm(method2 ~ method1)
  tmp.mean <- mean(tmp.ccc$blalt$delta)
  tmp.sd <- sqrt(var(tmp.ccc$blalt$delta))
  if (is.null(xrange))
    xrange <- range(method1)
  if (is.null(yrange))
    yrange <- range(method2)
  if (Ptype == "scatter") {  # Scatter Plot
    plot(method1, method2, xlab = xlabel, xlim = xrange, ylim = yrange,
         ylab = ylabel, pch = 16, sub = paste("(", subtitle, ")"))
    abline(a = 0, b = 1, lty = 2)
    abline(z, lty = 1)
    usr <- par("usr")  	# get user coordinates
    par(usr = c(0, 1, 0, 1))  # new relative user coordinates
    text(0.5, 0.18, r2lab, adj = 0)
    text(0.5, 0.23, r2rob, adj = 0)
    text(0.5, 0.12, Acc, adj = 0)
    text(0.5, 0.05, cclab, adj = 0)
    par(usr = usr)	# restore original user coordinates
  } else if (Ptype == "MAplot") {  # Bland-Altman or MAplot
    plot(tmp.ccc$blalt$mean, tmp.ccc$blalt$delta, pch = 16, xlab = "Average",
         ylab = "Difference", sub = paste("(", subtitle, ")"), ylim = MArange)
    abline(h = tmp.mean, lty = 1, col = "gray")
    abline(h = tmp.mean - (2 * tmp.sd), lty = 2, col = "gray")
    abline(h = tmp.mean + (2 * tmp.sd), lty = 2, col = "gray")
    abline(h = 0, lty = 1, col = "red")
    usr <- par("usr")  	# get user coordinates
    par(usr = c(0, 1, 0, 1))  # new relative user coordinates
    text(0.5, 0.95, r2lab, adj = 0)
    text(0.5, 0.88, Acc, adj = 0)
    text(0.5, 0.81, cclab, adj = 0)
    par(usr = usr)	# restore original user coordinates
  }
  if (metrics == TRUE) {
    return(c(Rc = round(tmp.ccc$rho.c[, 1], digits = 2),
             Ca = round(tmp.ccc$C.b, 2),
             R2 = round(cor(method1, method2), 2)))
  }
}
