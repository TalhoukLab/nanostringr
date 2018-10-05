#' Read NanoString RCC files
#'
#' Read a directory of RCC files and merge into count data and sample
#' attributes.
#'
#' RCC files for a sample are direct outputs from NanoString runs. We can
#' extract counts for each gene in a sample. Sample attributes include
#' information about Fov count, Fov counted, binding density, date, GeneRLF,
#' SampleID, etc. Both count and attribute data are then merged across samples.
#'
#' @param path directory path where RCC files exist
#' @param pattern file pattern for matching file extensions ".RCC" (or ".rcc")
#'
#' @return A tibble with columns "File.Name", "fov.count", "fov.counted",
#'   "binding.density".
#' @export
read_rcc <- function(path = ".", pattern = "\\.RCC$") {
  if (!dir.exists(path)) {
    utils::unzip(zipfile = paste0(path, ".ZIP"), exdir = path)
  }
  rcc_files <- list.files(
    path = path,
    pattern = pattern,
    full.names = TRUE,
    ignore.case = TRUE
  )
  raw <- rcc_files %>%
    purrr::map(parse_counts) %>%
    purrr::reduce(dplyr::inner_join, by = c("Code.Class", "Name", "Accession"))
  exp <- rcc_files %>%
    purrr::map_df(parse_attributes) %>%
    as.data.frame() %>%
    dplyr::mutate_at(c("fov.count", "fov.counted", "binding.density"),
                     as.numeric) %>%
    dplyr::mutate_at("nanostring.date",
                     dplyr::funs(as.character(as.Date(., "%Y%m%d"))))
  tibble::lst(raw, exp)
}

#' @param file RCC file path
#' @noRd
parse_counts <- function(file) {
  rcc_file <- readr::read_lines(file)
  sample_name <- get_attr(rcc_file, "^ID,.*[[:alpha:]]")
  cs_header <- grep("<Code_Summary>", rcc_file) + 1
  cs_last <- grep("</Code_Summary>", rcc_file) - 1
  rcc_parsed <-
    rcc_file[purrr::invoke(seq, c(cs_header, cs_last))] %>%
    paste(collapse = "\n") %>%
    readr::read_csv() %>%
    dplyr::rename(Code.Class = .data$CodeClass,
                  !!sample_name := .data$Count) %>%
    as.data.frame()
}

#' @param file RCC file path
#' @noRd
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
