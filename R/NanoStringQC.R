#' QC metrics for NanoString Data
#'
#' Computes and returns NanoString quality control metrics and flags.
#'
#' @param raw data frame of raw counts obtained from nCounter (rows represent
#'   genes, columns represent samples). The first three columns must be labeled:
#'   `c("Code.Class", "Name", "Accession")` and contain that information.
#' @param exp data frame of annotations with rows in the same order as the
#'   columns of `raw`. Requires a column labeled `"File.Name"` with entries
#'   corresponding to sample names in `raw`, also needs columns
#'   `c("fov.counted", "fov.count", "binding.density")`.These fields can be
#'   extracted from the nanostring RCC files.
#' @param detect threshold of percentage of genes expressed over limit of
#'   detection (LOD) that we would like to detect (not decimal), defaults to 80
#'   percent.
#' @param sn signal to noise ratio of the housekeeping genes we are willing to
#'   tolerate, defaults to 150.
#' @return data frame of annotations updated with normalization parameters
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

  # Extract PC gene concentrations
  PCgenes <- raw[raw$Code.Class == "Positive", "Name", drop = TRUE]
  if (!all(grepl("[[:digit:]]", PCgenes))) {
    stop("Positive controls need numeric concentrations: e.g. POS_A(128)")
  }
  PCconc <- as.numeric(regmatches(PCgenes, regexpr("\\d+\\.*\\d*", PCgenes)))

  # Expressions by code class
  endo_exp <- raw[raw$Code.Class == "Endogenous", -1:-3, drop = FALSE]
  neg_exp <- raw[raw$Code.Class == "Negative", -1:-3, drop = FALSE]
  pos_exp <- raw[raw$Code.Class == "Positive", -1:-3, drop = FALSE]
  hk_exp <- raw[raw$Code.Class == "Housekeeping", -1:-3, drop = FALSE]
  low_pos_exp <- as.list(raw[raw$Name %in% c("POS_E(0.5)", "POS_1"), -1:-3])

  # Code QC measures and flags
  exp |>
    dplyr::mutate(
      linPC = purrr::map_dbl(pos_exp, ~ round(summary(lm(. ~ PCconc))$r.squared, 2)),
      linFlag = .data$linPC < 0.95 | is.na(.data$linPC),
      perFOV = (.data$fov.counted / .data$fov.count) * 100,
      imagingFlag = .data$perFOV < 75,
      ncgMean = colMeans(neg_exp),
      ncgSD = purrr::map_dbl(neg_exp, sd),
      lod = .data$ncgMean + 2 * .data$ncgSD,
      llod = .data$ncgMean - 2 * .data$ncgSD,
      spcFlag = low_pos_exp < .data$llod | .data$ncgMean == 0,
      gd = purrr::map2_int(endo_exp, .data$lod, ~ sum(.x > .y)),
      pergd = (.data$gd / sum(raw$Code.Class == "Endogenous")) * 100,
      averageHK = purrr::map_dbl(hk_exp, ~ exp(mean(log2(.)))),
      snr = ifelse(.data$lod < 0.001, 0, .data$averageHK / .data$lod),
      bdFlag = .data$binding.density < 0.05 | .data$binding.density > 2.25,
      normFlag = .data$snr < sn | .data$pergd < detect,
      QCFlag = .data$linFlag | .data$imagingFlag | .data$spcFlag | .data$normFlag,
      dplyr::across(dplyr::matches("Flag"), ~ factor(
        .,
        levels = c(TRUE, FALSE),
        labels = c("Failed", "Passed")
      ))
    ) |>
    dplyr::select(-c("ncgMean", "ncgSD"))
}
