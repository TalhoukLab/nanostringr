#' QC for NanoString Data
#'
#' Computes NanoString quality control parameters and flags
#' @param raw matrix of raw counts obtained from nCounter (rows are genes).
#' The first three columns must be labeled: \code{c("Code.Class", "Name", "Accession")} and contain that information.
#' @param exp matrix of annotations with rows in the same order as the columns of raw. Needs a column labeled File.Name with entries corresponding to sample names in raw count, also needs columns fov.counted and fov.count as well as binding.density. These fields can be extracted from the RCC files.
#' @param detect percent threshold of genes over load that we would like to detect (not decimal).
#' @param sn Signal to noise ratio of the housekeeping genes we are willing to tolerate
#' @param plots logical; indicates whether plots to visualise the results are requested, defaults to \code{TRUE}
#' @param ttl a string to show the title on the plots
#' @param explore returns the plots only, defaults to \code{TRUE}
#' @return matrix of annotations updated with normalization parameters
#' @author Aline Talhouk, Derek Chiu
#' @import dplyr
#' @export
#' @examples
#' # Load otta package for raw datasets and annotation matrix
#' library(otta)
#' data(rawOVCA2, rawPROT, rawOTTA, annot)
#'
#' # Codeset 1, 2, 3 and annotations
#' cs1 <- rawOVCA2; cs2 <- rawPROT; cs3 <- rawOTTA; exp0 <- annot
#' exp0$geneRLF <- as.character(factor(exp0$geneRLF,
#'                      labels = c("HL1", "HL2", "HL3", "HuRef", "CS3", "mini", "CS1", "CS2")))
#'
#' # Compute NanoString QC
#' exp.CS1 <- NanoStringQC(cs1, exp0[exp0$geneRLF == "CS1", ],
#'                         plots = FALSE, detect = 50, ttl = "CodeSet 1")
#' exp.CS2 <- NanoStringQC(cs2, exp0[exp0$geneRLF == "CS2", ],
#'                         plots = FALSE, sn = 100, ttl = "CodeSet 2")
#' exp.CS3 <- NanoStringQC(cs3, exp0[exp0$geneRLF == "CS3", ],
#'                         plots = TRUE, detect = 50, sn = 100, ttl = "CodeSet 3")
NanoStringQC <- function(raw, exp, detect = 80, sn = 150, plots = TRUE,
                          ttl = " ", explore = TRUE) {
  
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

  if (plots) {
    par(mfrow = c(1, 2))
    plot(exp$sn, exp$pergd, pch = 20, col = "deepskyblue", xaxt = "n",
         ylim = c(0, 100), xlab = "Signal/Noise Ratio",
         ylab = "Percent of Genes Detected")
    axis(1, at = seq(0, max(exp$sn) + 1, 300))
    abline(v = sn, col = "red", lwd = 2)
    abline(h = detect, lty = 2)
    hist(exp$sn, 50, col = "cornsilk", prob = TRUE,
         xlab = "Signal to Noise", ylab = "Probability", main = "")
    abline(v = sn, lwd = 3)
    title(ttl, outer = TRUE, line = -2)
  }
  if (!explore)
    return(exp$QCFlag)
  else
    return(exp)
}
