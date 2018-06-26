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
    as.data.frame()
  tibble::lst(raw, exp)
}

#' @param file RCC file path
#' @noRd
parse_counts <- function(file) {
  rcc_file <- readr::read_lines(file)
  sample_name <- get_sample_name(rcc_file)
  cs_header <- grep("<Code_Summary>", rcc_file) + 1
  cs_last <- grep("</Code_Summary>", rcc_file) - 1
  rcc_parsed <-
    rcc_file[purrr::invoke(seq, c(cs_header, cs_last))] %>%
    readr::read_csv() %>%
    dplyr::rename(Code.Class = .data$CodeClass,
                  !!sample_name := .data$Count) %>%
    as.data.frame()
}

#' @param file RCC file path
#' @noRd
parse_attributes <- function(file) {
  rcc_file <- readr::read_lines(file)
  File.Name <- get_sample_name(rcc_file)
  fov.count <- get_attr(rcc_file, "FovCount")
  fov.counted <- get_attr(rcc_file, "FovCounted")
  binding.density <- get_attr(rcc_file, "BindingDensity")
  tibble::lst(File.Name, fov.count, fov.counted, binding.density)
}

#' @param rcc_file RCC file
#' @param attr name of lane attribute to extract
#' @noRd
get_attr <- function(rcc_file, attr) {
  grep(paste0(attr, ","), rcc_file, value = TRUE) %>%
    strsplit(split = ",") %>%
    purrr::pluck(1, 2) %>%
    as.numeric()
}

#' @param rcc_file RCC file
#' @noRd
get_sample_name <- function(rcc_file) {
  id <- grep("<Sample_Attributes>", rcc_file) + 1
  rcc_file[id] %>%
    strsplit(split = ",") %>%
    purrr::pluck(1, 2) %>%
    make.names()
}
