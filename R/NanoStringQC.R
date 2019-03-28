#' QC metrics for NanoString Data
#'
#' Computes and returns NanoString quality control metrics and flags.
#'
#' @param raw matrix of raw counts obtained from nCounter (rows represent genes,
#'   columns represent samples). The first three columns must be labeled:
#'   `c("Code.Class", "Name", "Accession")` and contain that information.
#' @param exp matrix of annotations with rows in the same order as the columns
#'   of `raw`. Requires a column labeled `"File.Name"` with entries
#'   corresponding to sample names in `raw`, also needs columns
#'   `c("fov.counted", "fov.count", "binding.density")`.These fields can be
#'   extracted from the nanostring RCC files.
#' @param detect threshold of percentage of genes expressed over limit of
#'   detection (LOD) that we would like to detect (not decimal), defaults to 80
#'   percent.
#' @param sn signal to noise ratio of the housekeeping genes we are willing to
#'   tolerate, defaults to 150.
#' @return matrix of annotations updated with normalization parameters
#' @author Aline Talhouk, Derek Chiu
#' @export
#' @examples
#' exp.OVD <- subset(expQC, OVD == "Yes")
#' expOVD <- NanoStringQC(ovd.r, exp.OVD)
NanoStringQC <- function(raw, exp, detect = 80, sn = 150) {
  # Run checks to make sure the data is in the right order
  assertthat::assert_that(check_colnames(raw))  # Checks format of raw counts
  assertthat::assert_that(check_genes(raw))  # Checks HK genes are specified
  assertthat::assert_that(ncol(raw) == nrow(exp) + 3)  # Checks data dimensions

  sn.in <- sn
  genes <- raw$Name
  rownames(raw) <- genes
  HKgenes <- genes[raw$Code.Class == "Housekeeping"]
  PCgenes <- genes[raw$Code.Class == "Positive"]
  NCgenes <- genes[raw$Code.Class == "Negative"]
  Hybgenes <- genes[raw$Code.Class != "Endogenous"]
  if (!all(grepl("[[:digit:]]", PCgenes))) {
    stop("Positive controls must have concentrations in brackets: ex POS_A(128)")
  }
  PCconc <- as.numeric(gsub(".*\\((.*)\\).*", "\\1", PCgenes))
  flag.levs <- c("Failed", "Passed")
  exp %>%
    dplyr::mutate(
      linPC = raw[PCgenes, -(1:3), drop = FALSE] %>%
        purrr::map_dbl(~ summary(lm(. ~ PCconc))$r.squared) %>%
        round(2),
      linFlag = factor(ifelse(.data$linPC < 0.95 | is.na(.data$linPC), "Failed", "Passed"),
                       flag.levs),
      perFOV = (.data$fov.counted / .data$fov.count) * 100,
      imagingFlag = factor(ifelse(.data$perFOV < 75, "Failed", "Passed"), flag.levs),
      ncgMean = purrr::map_dbl(raw[NCgenes, -(1:3), drop = FALSE], mean),
      ncgSD = purrr::map_dbl(raw[NCgenes, -(1:3), drop = FALSE], sd),
      lod = .data$ncgMean + 2 * .data$ncgSD,
      llod = .data$ncgMean - 2 * .data$ncgSD,
      spcFlag = factor(ifelse(
        t(as.vector(raw["POS_E(0.5)", -(1:3), drop = FALSE]) < .data$llod |
            .data$ncgMean == 0),
        "Failed", "Passed"), flag.levs),
      gd = apply(raw[!(rownames(raw) %in% Hybgenes), -(1:3), drop = FALSE] > .data$lod, 2, sum),
      pergd = (.data$gd / nrow(raw[!(rownames(raw) %in% Hybgenes), -(1:3), drop = FALSE])) * 100,
      averageHK = exp(purrr::map_dbl(log2(raw[HKgenes, -(1:3), drop = FALSE]), mean)),
      sn = ifelse(.data$lod < 0.001, 0, .data$averageHK / .data$lod),
      bdFlag = factor(ifelse(
        .data$binding.density < 0.05 | .data$binding.density > 2.25,
        "Failed", "Passed"), flag.levs),
      normFlag = factor(ifelse(
        .data$sn < sn.in | .data$pergd < detect,
        "Failed", "Passed"), flag.levs),
      QCFlag = factor(ifelse(
        .data$spcFlag == "Failed" | .data$imagingFlag == "Failed" | .data$linFlag == "Failed",
        "Failed", "Passed"), flag.levs)
    ) %>%
    magrittr::set_rownames(rownames(exp)) %>%
    dplyr::select(-c(.data$ncgMean, .data$ncgSD))
}
