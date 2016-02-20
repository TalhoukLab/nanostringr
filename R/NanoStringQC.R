#' QC metrics for NanoString Data
#'
#' Computes and returns NanoString quality control metrics and flags.
#'
#' @param raw matrix of raw counts obtained from nCounter (rows represent genes, columns represent samples).
#' The first three columns must be labeled: \code{c("Code.Class", "Name", "Accession")} and contain that information.
#' @param exp matrix of annotations with rows in the same order as the columns of \code{raw}. Requires a column labeled \code{"File.Name"} with entries corresponding to sample names in \code{raw}, also needs columns \code{c("fov.counted", "fov.count", "binding.density")}.These fields can be extracted from the nanostring RCC files.
#' @param detect threshold of percentage of genes expressed over limit of detection (LOD) that we would like to detect (not decimal), defaults to 80 percent.
#' @param sn signal to noise ratio of the housekeeping genes we are willing to tolerate, defaults to 150.
#'
#' @return matrix of annotations updated with normalization parameters
#'
#' @author Aline Talhouk, Derek Chiu
#'
#' @import dplyr
#'
#' @export
#'
#' @examples
#' exp.OVD <-subset(expQC,OVD=="Yes")
#' expOVD <- NanoStringQC(ovd.r,exp.OVD)


NanoStringQC <- function(raw, exp, detect = 80, sn = 150) {

  # Run a bunch of checks to make sure the data is in the right order
  assertthat::assert_that(check_colnames(raw)) #Checks format of raw counts
  assertthat::assert_that(check_genes(raw)) #Checks that HK genes are specified
  assertthat::assert_that(ncol(raw) == nrow(exp) + 3)
  assertthat::assert_that(all(substring(colnames(raw[,-(1:3)]),2) == exp$File.Name))

  sn.in <- sn
  genes <- raw$Name
  rownames(raw) <- genes
  HKgenes <- genes[raw$Code.Class == "Housekeeping"]
  PCgenes <- genes[raw$Code.Class == "Positive"]
  NCgenes <- genes[raw$Code.Class == "Negative"]
  Hybgenes <- genes[raw$Code.Class != "Endogenous"]
  PCconc <- as.numeric(sub("\\).*", "", sub(".*\\(", "", PCgenes)))
  flag.levs <- c("Failed", "Passed")
  . <- linPC <- linFlag <- fov.counted <- fov.count <- perFOV <- ncgMean <-
    ncgSD <- llod <- lod <- gd <- averageHK <- binding.density <- pergd <-
    spcFlag <- normFlag <- imagingFlag <- linFlag <- rn <- NULL

  exp <- exp %>%
    mutate(rn = rownames(.),
           linPC = round(apply(raw[PCgenes, -(1:3)], 2,
                               function(x)
                                 summary(lm(x ~ PCconc))$r.squared), 2),
           linFlag = factor(ifelse(linPC < 0.95 | is.na(linPC),
                                   "Failed", "Passed"), flag.levs),
           perFOV = (fov.counted / fov.count) * 100,
           imagingFlag = factor(ifelse(perFOV < 75, "Failed", "Passed"),
                                flag.levs),
           ncgMean = apply(raw[NCgenes, -(1:3)], 2, mean),
           ncgSD = apply(raw[NCgenes, -(1:3)], 2, sd),
           lod = ncgMean + 2 * ncgSD,
           llod = ncgMean - 2 * ncgSD,
           spcFlag = factor(
             ifelse(t(as.vector(raw["POS_E(0.5)", -(1:3)]) < llod |
                        ncgMean == 0), "Failed", "Passed"), flag.levs),
           gd = apply(raw[!(rownames(raw) %in% Hybgenes), -(1:3)] > lod,
                      2, sum),
           pergd = (gd / nrow(raw[!(rownames(raw) %in% Hybgenes),
                                  -(1:3)])) * 100,
           averageHK = exp(apply(log2(raw[HKgenes, -(1:3)]), 2, mean)),
           sn = ifelse(lod < 0.001, 0, averageHK / lod),
           bdFlag = factor(
             ifelse(binding.density < 0.05 | binding.density > 2.25,
                    "Failed", "Passed"), flag.levs),
           normFlag = factor(
             ifelse(sn < sn.in | pergd < detect,
                    "Failed", "Passed"), flag.levs),
           QCFlag = factor(ifelse(
             as.vector(spcFlag == "Failed"|
                         imagingFlag == "Failed" | linFlag == "Failed"),
             "Failed", "Passed"))) %>%
    magrittr::set_rownames(.$rn) %>%
    select(-rn, -ncgMean, -ncgSD)
  names(exp[, ]) <- NULL
  return(exp)
}
