# magrittr placeholder
globalVariables(".")

#' Join grouping dataset to CodeSet
#'
#' Join grouping dataset to CodeSet, then keep or discard the common ids, and
#' finally average the duplicate samples. Used as a helper function for
#' `normalize_random()`.
#'
#' @param x nanostring expression dataset
#' @param y grouping dataset (e.g. random samples from each histotype)
#' @param id identifier to join `x` and `y` by
#' @param type the type of join to perform. If "keep", then a semi join keeps
#'   cases in `x` that also appear in `y`. If "discard", then an anti join
#'   discards cases in `x` that appear in `y`.
#' @author Derek Chiu
#' @noRd
join_avg <- function(x, y, id, type = c("keep", "discard")) {
  type <- match.arg(type)
  join_fun <- switch(type,
                     keep = dplyr::semi_join,
                     discard = dplyr::anti_join)
  x |>
    join_fun(y, by = id) |>
    dplyr::summarize(dplyr::across(dplyr::where(is.double), mean),
                     .by = !!rlang::sym(id)) |>
    dplyr::arrange(!!rlang::sym(id)) |>
    tibble::column_to_rownames("ottaID")
}

#' Concordance correlation coefficient
#'
#' Using [epiR::epi.ccc()] to avoid installing downstream dependencies of
#' `epiR` like `sf` that may be difficult to install from source.
#'
#' @noRd
epi_ccc <- function (x, y, ci = "z-transform", conf.level = 0.95, rep.measure = FALSE,
                     subjectid) {
  N. <- 1 - ((1 - conf.level)/2)
  zv <- stats::qnorm(N., mean = 0, sd = 1)
  dat <- data.frame(x, y)
  id <- stats::complete.cases(dat)
  nmissing <- sum(!stats::complete.cases(dat))
  dat <- dat[id, ]
  k <- length(dat$y)
  yb <- mean(dat$y)
  sy2 <- var(dat$y) * (k - 1)/k
  sd1 <- sd(dat$y)
  xb <- mean(dat$x)
  sx2 <- var(dat$x) * (k - 1)/k
  sd2 <- sd(dat$x)
  r <- cor(dat$x, dat$y)
  sl <- r * sd1/sd2
  sxy <- r * sqrt(sx2 * sy2)
  p <- 2 * sxy/(sx2 + sy2 + (yb - xb)^2)
  delta <- (dat$x - dat$y)
  rmean <- apply(dat, MARGIN = 1, FUN = mean)
  blalt <- data.frame(mean = rmean, delta)
  v <- sd1/sd2
  u <- (yb - xb)/((sx2 * sy2)^0.25)
  C.b <- p/r
  sep <- sqrt(((1 - ((r)^2)) * (p)^2 * (1 - ((p)^2))/(r)^2 +
                 (2 * (p)^3 * (1 - p) * (u)^2/r) - 0.5 * (p)^4 * (u)^4/(r)^2)/(k - 2))
  ll <- p - (zv * sep)
  ul <- p + (zv * sep)
  t <- log((1 + p)/(1 - p))/2
  set = sep/(1 - ((p)^2))
  llt = t - (zv * set)
  ult = t + (zv * set)
  llt = (exp(2 * llt) - 1)/(exp(2 * llt) + 1)
  ult = (exp(2 * ult) - 1)/(exp(2 * ult) + 1)
  if (rep.measure == TRUE) {
    dat$sub <- subjectid
    if (!is.factor(dat$sub))
      dat$sub <- as.factor(dat$sub)
    nsub <- length(levels(dat$sub))
    model <- stats::aov(delta ~ dat$sub)
    MSB <- stats::anova(model)[[3]][1]
    MSW <- stats::anova(model)[[3]][2]
    pairs <- NULL
    for (i in 1:nsub) {
      pairs[i] <- sum(is.na(delta[dat$sub == levels(dat$sub)[i]]) ==
                        FALSE)
    }
    sig.dl <- (MSB - MSW)/((sum(pairs)^2 - sum(pairs^2))/((nsub -
                                                             1) * sum(pairs)))
    delta.sd <- sqrt(sig.dl + MSW)
  }
  if (rep.measure == FALSE) {
    delta.sd <- sqrt(var(delta, na.rm = TRUE))
  }
  ba.p <- mean(delta)
  ba.l <- ba.p - (zv * delta.sd)
  ba.u <- ba.p + (zv * delta.sd)
  sblalt <- data.frame(est = ba.p, delta.sd = delta.sd, lower = ba.l,
                       upper = ba.u)
  if (ci == "asymptotic") {
    rho.c <- data.frame(p, ll, ul)
    names(rho.c) <- c("est", "lower", "upper")
    rval <- list(rho.c = rho.c, s.shift = v, l.shift = u,
                 C.b = C.b, blalt = blalt, sblalt = sblalt, nmissing = nmissing)
  }
  else if (ci == "z-transform") {
    rho.c <- data.frame(p, llt, ult)
    names(rho.c) <- c("est", "lower", "upper")
    rval <- list(rho.c = rho.c, s.shift = v, l.shift = u,
                 C.b = C.b, blalt = blalt, sblalt = sblalt, nmissing = nmissing)
  }
  return(rval)
}
