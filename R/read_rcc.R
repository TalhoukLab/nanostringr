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
#' read_rcc(dirname(rcc_file))
read_rcc <- function(path = ".") {
  if (!dir.exists(path)) {
    utils::unzip(zipfile = paste0(path, ".ZIP"), exdir = path)
  }
  rcc_files <- list.files(path,
                          pattern = "\\.RCC$",
                          full.names = TRUE,
                          ignore.case = TRUE)
  rcc_basenames <- tools::file_path_sans_ext(basename(rcc_files))
  rcc_files <- stats::setNames(rcc_files, rcc_basenames)
  raw <- rcc_files |>
    purrr::map(parse_counts) |>
    purrr::imap(~ dplyr::rename_with(.x, rlang::quo(.y), 4)) |>
    purrr::reduce(\(x, y) dplyr::inner_join(x, y, by = c("Code.Class", "Name", "Accession"))) |>
    dplyr::mutate(!!"Name" := dplyr::case_match(
      .data$Name,
      "CD3E" ~ "CD3e",
      "PD1" ~ "PD-1",
      "PDL1" ~ "PD-L1",
      .default = .data$Name
    )) |>
    as.data.frame()
  exp <- rcc_files |>
    purrr::map(parse_attributes) |>
    dplyr::bind_rows(.id = "sample") |>
    dplyr::select(-"File.Name") |>
    as.data.frame()
  rlang::dots_list(raw, exp, .named = TRUE)
}

#' @param file RCC file name
#' @name rcc
#' @return `parse_counts()` reads a single RCC file and returns a tibble of
#'   parsed counts.
#' @export
parse_counts <- function(file) {
  rcc_file <- readLines(file)
  sample_name <- get_attr(rcc_file, "^ID,.*[[:alnum:]]")
  cs_header <- grep("<Code_Summary>", rcc_file) + 1
  cs_last <- grep("</Code_Summary>", rcc_file) - 1
  rcc_file[cs_header:cs_last] |>
    paste(collapse = "\n") |>
    utils::read.csv(text = _) |>
    dplyr::rename(Code.Class = "CodeClass", !!sample_name := "Count") |>
    dplyr::as_tibble()
}

#' @name rcc
#' @return `parse_attributes()` reads a single RCC file and returns a list of
#'   parsed attributes.
#' @export
parse_attributes <- function(file) {
  rcc_file <- readLines(file)
  rcc_attrs <- c(
    File.Name = "^ID,.*[[:alnum:]]",
    geneRLF = "GeneRLF",
    nanostring.date = "Date",
    cartridgeID = "CartridgeID",
    lane.number = "^ID,[[:digit:]]+$",
    fov.count = "FovCount",
    fov.counted = "FovCounted",
    binding.density =  "BindingDensity"
  )
  rcc_attrs |>
    purrr::map(get_attr, rcc_file = rcc_file) |>
    purrr::map_at(c("fov.count", "fov.counted", "binding.density"),
                  as.numeric) |>
    purrr::map_at("nanostring.date", ~ as.character(as.Date(., "%Y%m%d")))
}

#' @param rcc_file RCC file
#' @param attr name of lane attribute to extract
#' @noRd
get_attr <- function(rcc_file, attr) {
  grep(attr, rcc_file, value = TRUE) |>
    strsplit(split = ",") |>
    purrr::pluck(1, 2)
}
