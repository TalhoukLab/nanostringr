#' Read NanoString RCC files
#'
#' Read RCC files and extract count and attribute data. Use `read_rcc()` for
#' multiple files, and use the `parse_*()` functions for single files.
#'
#' RCC files for a sample are direct outputs from NanoString runs. We can
#' extract counts for each gene in a sample. Sample attributes include sample
#' ID, GeneRLF, date, cartridge ID, lane number, Fov count, Fov counted, and
#' binding density. `read_rcc()` merges both count and attribute data across
#' samples.
#'
#' If `path` points to a zipped RCC file with multiple samples, the zip file is
#' uncompressed and a directory of RCC sample files is created with the same
#' name. Only file extensions ".RCC" or ".rcc" are allowed.
#'
#' @param path directory path for multiple RCC files
#'
#' @return `read_rcc()` reads in a directory of RCC files and outputs a list
#'   with two elements:
#' * `raw`: A tibble of parsed counts for multiple RCC files created by calling
#' `parse_counts()` on each sample. Columns include "Code.Class", "Name",
#' "Accession", and a column for each sample ID. There is one row per gene.
#' * `exp`: A tibble of parsed attributes for multiple RCC files created by
#' calling `parse_attributes()` on each sample. Columns include "File.Name"
#' (sample ID), "geneRLF", "nanostring.date", "cartridgeID", "lane.number",
#' fov.count", "fov.counted", "binding.density". There is one row per sample.
#'
#' @author Derek Chiu
#' @name rcc
#' @export
#' @examples
#' rcc_file <- system.file("extdata", "example.RCC", package = "nanostringr")
#' parse_counts(rcc_file)
#' parse_attributes(rcc_file)
read_rcc <- function(path = ".") {
  if (!dir.exists(path)) {
    utils::unzip(zipfile = paste0(path, ".ZIP"), exdir = path)
  }
  rcc_files <-
    list.files(path, pattern = "\\.RCC$", full.names = TRUE, ignore.case = TRUE)
  raw <- rcc_files %>%
    purrr::map(parse_counts) %>%
    purrr::reduce(dplyr::inner_join, by = c("Code.Class", "Name", "Accession"))
  exp <- rcc_files %>%
    purrr::map_df(parse_attributes)
  tibble::lst(raw, exp)
}

#' @param file RCC file name
#' @name rcc
#' @return `parse_counts()` reads a single RCC file and returns a tibble of
#'   parsed counts.
#' @export
parse_counts <- function(file) {
  rcc_file <- readr::read_lines(file)
  sample_name <- get_attr(rcc_file, "^ID,.*[[:alnum:]]")
  cs_header <- grep("<Code_Summary>", rcc_file) + 1
  cs_last <- grep("</Code_Summary>", rcc_file) - 1
  rcc_file[cs_header:cs_last] %>%
    paste(collapse = "\n") %>%
    readr::read_csv() %>%
    dplyr::rename(Code.Class = .data$CodeClass, !!sample_name := .data$Count)
}

#' @inheritParams parse_counts
#' @name rcc
#' @return `parse_attributes()` reads a single RCC file and returns a list of
#'   parsed attributes.
#' @export
parse_attributes <- function(file) {
  rcc_file <- readr::read_lines(file)
  attr_patterns <- c("^ID,.*[[:alnum:]]", "GeneRLF", "Date", "CartridgeID",
                     "^ID,[[:digit:]]+$", "FovCount", "FovCounted",
                     "BindingDensity")
  attr_names <- c("File.Name", "geneRLF", "nanostring.date", "cartridgeID",
                  "lane.number", "fov.count", "fov.counted", "binding.density")
  attr_patterns %>%
    purrr::set_names(attr_names) %>%
    purrr::map(get_attr, rcc_file = rcc_file) %>%
    purrr::map_at(c("fov.count", "fov.counted", "binding.density"),
                  as.numeric) %>%
    purrr::map_at("nanostring.date", ~ as.character(as.Date(., "%Y%m%d")))
}

#' @param rcc_file RCC file
#' @param attr name of lane attribute to extract
#' @noRd
get_attr <- function(rcc_file, attr) {
  grep(attr, rcc_file, value = TRUE) %>%
    strsplit(split = ",") %>%
    purrr::pluck(1, 2)
}
