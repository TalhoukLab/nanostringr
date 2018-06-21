#' Convert directory of RCC files to an expression object
#'
#' Reads in RCC files from a directory and reformats to an expression object
#' used for analysis
#'
#' @param dir directory path where RCC files exist
#'
#' @return A tibble with columns "File.Name", "fov.count", "fov.counted",
#'   "binding.density".
#' @export
#' @examples
#' rcc_dir <- "~/Downloads/20180524_CAP run3 24May2018_RCC"
#' raw <- rcc_to_codeset(rcc_dir)
#' exp <- rcc_to_exp(rcc_dir)
#' NanoStringQC(raw, exp)
rcc_to_exp <- function(dir) {
  rcc_files <- list.files(dir, pattern = "RCC", full.names = TRUE)
  rcc_parsed_list <- purrr::map(rcc_files, parse_exp)
  as.data.frame(dplyr::bind_rows(rcc_parsed_list))
}

#' @param file RCC file path
#' @noRd
parse_exp <- function(file) {
  rcc_file <- readr::read_lines(file)
  File.id <- grep("<Sample_Attributes>", rcc_file) + 1
  File.Name <- make.names(strsplit(rcc_file[sample_id], ",")[[1]][2])
  fov.count <- get_attr(rcc_file, "FovCount")
  fov.counted <- get_attr(rcc_file, "FovCounted")
  binding.density <- get_attr(rcc_file, "BindingDensity")
  tibble::lst(File.Name, fov.count, fov.counted, binding.density)
}

#' @param attr name of lane attribute to extract
#' @noRd
get_attr <- function(rcc_file, attr) {
  grep(paste0(attr, ","), rcc_file, value = TRUE) %>%
    strsplit(split = ",") %>%
    purrr::pluck(1, 2) %>%
    as.numeric()
}
