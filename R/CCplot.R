#' Concordance Correlation Plot
#'
#' Plotting function for reliability measure.
#'
#' @param method1 measurements obtained in batch 1 or using method 1
#' @param method2 measurements obtained in batch 2 or using method 2
#' @param Ptype type of plot to be outputted c("scatter", "MAplot")
#' @param metrics if `TRUE`, prints Rc, Ca, and R2 to console
#' @param xlabel x-axis label for scatterplot
#' @param ylabel y-axis label for scatterplot
#' @param title title for the main plot
#' @param subtitle subtitle of plot
#' @param xrange range of x axis
#' @param yrange range of y axis
#' @param MArange MA range
#' @return Either a scatterplot or MA plot showing concordance correlation.
#' @author Aline Talhouk
#' @importFrom graphics abline par plot text
#' @importFrom stats cor lm sd var
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
#'
#' # One scatterplot
#' CCplot(a1, a2, Ptype = "scatter")
#'
#'
#' m2 <- list(a1, a2, a3, a4, a5, a6)
#' mains <- c("Perfect Agreement", "Very Good Agreement", "Location Shift",
#'            "Scale Shift", "Location and Scale Shift", "Measurement Error")
#' subs <- letters[1:6]
#' par(mfrow = c(3, 2), mar = c(5.1, 4.1, 1.5, 1.5))
#'
#' # Scatterplots
#' mapply(function(y, t, s)
#'   CCplot(method1 = a1, method2 = y, Ptype = "scatter",
#'          xlabel = "X", ylabel = "Y", title = t, subtitle = s),
#'   y = m2, t = mains, s = subs)
#'
#' # MAplots and show metrics
#' mapply(function(y, t, s)
#'   CCplot(method1 = a1, method2 = y, Ptype = "MAplot",
#'          title = t, subtitle = s, metrics = TRUE),
#'   y = m2, t = mains, s = subs)
CCplot <- function(method1, method2, Ptype = "None", metrics = FALSE,
                   xlabel = "", ylabel = "", title = "", subtitle = NULL,
                   xrange = NULL, yrange = NULL, MArange = c(-3.5, 5.5)) {
  assertthat::assert_that(length(method1) == length(method2))
  tmp.ccc <- epiR::epi.ccc(method1, method2, ci = "z-transform",
                           conf.level = 0.95)
  r2rob <- paste0("rM: ", round(ccaPP::corM(method1, method2), 2))
  r2lab <- paste0("R: ", round(cor(method1, method2), 2))
  Acc <- paste0("Ca: ", round(tmp.ccc$C.b, 2))
  cc <- round(tmp.ccc$rho.c, 2)
  cclab <- paste0("Rc: ", cc[1], " (", cc[2], " - ", cc[3], ")")
  z <- lm(method2 ~ method1)
  tmp.mean <- mean(tmp.ccc$blalt$delta)
  tmp.sd <- sd(tmp.ccc$blalt$delta)
  if (!is.null(subtitle)) {
    sub <- paste("(", subtitle, ")")
  } else {
    sub <- NULL
  }
  xrange <- xrange %||% range(method1)
  yrange <- yrange %||% range(method2)
  if (Ptype == "scatter") {  # Scatter Plot
    plot(method1, method2, xlab = xlabel, xlim = xrange, ylim = yrange,
         ylab = ylabel, main = title, pch = 16, sub = sub)
    abline(a = 0, b = 1, lty = 2)
    abline(z, lty = 1)
    usr <- par("usr")  	# get user coordinates
    par(usr = c(0, 1, 0, 1))  # new relative user coordinates
    text(0.5, 0.23, r2rob, adj = 0)
    text(0.5, 0.18, r2lab, adj = 0)
    text(0.5, 0.12, Acc, adj = 0)
    text(0.5, 0.05, cclab, adj = 0)
    par(usr = usr)	# restore original user coordinates
  } else if (Ptype == "MAplot") {  # Bland-Altman or MAplot
    plot(tmp.ccc$blalt$mean, tmp.ccc$blalt$delta, pch = 16, xlab = "Average",
         ylab = "Difference", main = title, sub = sub, ylim = MArange)
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
  if (metrics) {
    round(c(
      Rc = tmp.ccc$rho.c[, 1],
      Ca = tmp.ccc$C.b,
      R2 = cor(method1, method2)
    ), 2)
  }
}
