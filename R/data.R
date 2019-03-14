#' Expression QC data
#'
#' Quality control metrics for the five cohorts analyzed in NanoString
#' experiments.
#'
#' The total number of samples from the five cohorts is 561.
#'
#' @format A data frame with 561 rows and 23 columns.
#' @name expQC
#' @seealso [cohort]
NULL

#' NanoString Experiment Cohorts
#'
#' There were five different cohorts used in NanoString experiments.
#'
#' Each data object contains raw expression counts, so no normalization has been
#' applied. The format is a data frame with genes as rows, samples as columns.
#' Note that the first three columns contain gene metadata and are always
#' labelled "Code.Class", "Name", and "Accession", and the rest are sample
#' names. Hence, for the `hld.r` data, the raw counts are contained in 232 genes
#' for 77 - 3 = 74 samples. The total number of samples is 74 + 258 + 26 + 68 +
#' 135 = 561, which matches the number of rows in [expQC], the expression QC
#' data.
#'
#' @name cohort
#' @aliases hld.r ovd.r ovc.r hlo.r ovo.r
#' @format
#' * `hld.r` Hodgkin Lymphoma Clinical Samples: a data frame with 232 rows and
#' 77 columns
#' * `ovd.r` Ovarian Cancer Clinical Samples: a data frame with 133 rows and 261
#' columns
#' * `ovc.r` Ovarian Cancer Cell Lines: a data frame with 133 rows and 29
#' columns
#' * `hlo.r` DNA Oligonucleotides for the HL CodeSet: a data frame with 40 rows
#' and 71 columns
#' * `ovo.r` DNA Oligonucleotides for the OC CodeSet: a data frame with 133 rows
#' and 138 columns
#' @seealso [expQC]
#' @source See Table 1 of
#'   <https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0153844>
#'   for details.
NULL

#' @rdname cohort
"hld.r"

#' @rdname cohort
"ovd.r"

#' @rdname cohort
"ovc.r"

#' @rdname cohort
"hlo.r"

#' @rdname cohort
"ovo.r"
