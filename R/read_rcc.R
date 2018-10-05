#' Read NanoString RCC files
#'
#' Read directory of RCC files and merge samples into count and attribute
#' data.
#'
#' RCC files for a sample are direct outputs from NanoString runs. We can
#' extract counts for each gene in a sample. Sample attributes include sample
#' ID, GeneRLF, date, cartridge ID, lane number, Fov count, Fov counted, and
#' binding density. Both count and attribute data are then merged across
#' samples.
#'
#' If `path` points to a zipped RCC file with multiple samples, the zip file is
#' unzipped and a directory of RCC sample files are created with the same name.
#'
#' @param path directory path for multiple RCC files
#' @param pattern pattern for matching file extensions ".RCC" or ".rcc"
#'
#' @return A tibble with columns "File.Name", "fov.count", "fov.counted",
#'   "binding.density".
#' @author Derek Chiu
#' @name rcc
#' @export
read_rcc <- function(path = ".", pattern = "\\.RCC$") {
  if (!dir.exists(path)) {
    utils::unzip(zipfile = paste0(path, ".ZIP"), exdir = path)
  }
  rcc_files <-
    list.files(path, pattern, full.names = TRUE, ignore.case = TRUE)
  raw <- rcc_files %>%
    purrr::map(parse_counts) %>%
    purrr::reduce(dplyr::inner_join, by = c("Code.Class", "Name", "Accession"))
  exp <- rcc_files %>%
    purrr::map_df(parse_attributes) %>%
    dplyr::mutate_at(c("fov.count", "fov.counted", "binding.density"),
                     as.numeric) %>%
    dplyr::mutate(!!"nanostring.date" :=
                    as.character(as.Date(.data$nanostring.date, "%Y%m%d")))
  tibble::lst(raw, exp)
}

#' @param file RCC file name
#' @name rcc
#' @export
parse_counts <- function(file) {
  rcc_file <- readr::read_lines(file)
  sample_name <- get_attr(rcc_file, "^ID,.*[[:alpha:]]")
  cs_header <- grep("<Code_Summary>", rcc_file) + 1
  cs_last <- grep("</Code_Summary>", rcc_file) - 1
  rcc_parsed <-
    rcc_file[purrr::invoke(seq, c(cs_header, cs_last))] %>%
    paste(collapse = "\n") %>%
    readr::read_csv() %>%
    dplyr::rename(Code.Class = .data$CodeClass, !!sample_name := .data$Count)
}

#' @inheritParams parse_counts
#' @name rcc
#' @export
parse_attributes <- function(file) {
  rcc_file <- readr::read_lines(file)
  attr_patterns <- c("^ID,.*[[:alpha:]]", "GeneRLF", "Date", "CartridgeID",
                     "^ID,[[:digit:]]+$", "FovCount", "FovCounted",
                     "BindingDensity")
  attr_names <- c("File.Name", "geneRLF", "nanostring.date", "cartridgeID",
                  "lane.number", "fov.count", "fov.counted", "binding.density")
  attr_patterns %>%
    purrr::set_names(attr_names) %>%
    purrr::map(get_attr, rcc_file = rcc_file)
}

#' @param rcc_file RCC file
#' @param attr name of lane attribute to extract
#' @noRd
get_attr <- function(rcc_file, attr) {
  grep(attr, rcc_file, value = TRUE) %>%
    strsplit(split = ",") %>%
    purrr::pluck(1, 2)
}
