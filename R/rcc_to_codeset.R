#' Convert directory of RCC files to a codeset
#'
#' Reads in RCC files from a directory and reformats to a codeset that can be
#' used for analysis
#'
#' @param dir directory path where RCC files exist
#'
#' @return A tibble with columns "Code.Class", "Name", "Accession", and counts
#'   from each sample.
#' @export
rcc_to_codeset <- function(dir) {
  rcc_files <- list.files(dir, pattern = "\\.RCC$", full.names = TRUE)
  rcc_parsed_list <- purrr::map(rcc_files, parse_codeset)
  rcc_parsed_list %>%
    purrr::reduce(dplyr::inner_join,
                  by = c("Code.Class", "Name", "Accession")) %>%
    as.data.frame()
}

#' @param file RCC file path
#' @noRd
parse_codeset <- function(file) {
  rcc_file <- readr::read_lines(file)
  sample_id <- grep("<Sample_Attributes>", rcc_file) + 1
  sample_name <- make.names(strsplit(rcc_file[sample_id], ",")[[1]][2])
  cs_header <- grep("<Code_Summary>", rcc_file) + 1
  cs_last <- grep("</Code_Summary>", rcc_file) - 1
  rcc_parsed <-
    rcc_file[purrr::invoke(seq, c(cs_header, cs_last))] %>%
    readr::read_csv() %>%
    dplyr::rename(Code.Class = .data$CodeClass,
                  !!sample_name := .data$Count)
}
