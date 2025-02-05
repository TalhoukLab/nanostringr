nanostringr
===========

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/nanostringr)](https://cran.r-project.org/package=nanostringr)
[![R-CMD-check](https://github.com/TalhoukLab/nanostringr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/TalhoukLab/nanostringr/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/TalhoukLab/nanostringr/graph/badge.svg)](https://app.codecov.io/gh/TalhoukLab/nanostringr)
<!-- badges: end -->

Overview
--------

An R Package for quality assurance checking, normalization and batch effects adjustments of NanoString data suitable for single sample processing. This is the companion R package for the paper ["Single-Patient Molecular Testing with NanoString nCounter Data Using a Reference-Based Strategy for Batch Effect Correction"](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0153844), published in PLOS ONE.


Installation
------------

You can install `nanostringr` from CRAN with:

``` r
install.packages("nanostringr")
```

Or get the latest development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("TalhoukLab/nanostringr")
```


Usage
--------

To see the full list of exported functions:

``` r
library(nanostringr)
ls("package:nanostringr")
```

A quick overview of the key functions:

-   `NanoStringQC`: Computes quality assurance metrics.
-   `HKnorm`: Performs log (base 2) transformation and normalization to Housekeeping Genes
-   `refMethod`: Performs batch effect correction using the reference-based strategy
-   `parse_counts`: Read RCC files and extract raw counts data
-   `parse_attributes`: Read RCC files and extract expression annotation data

Please see `citation("nanostringr")` for information on how to cite the PLOS ONE paper when using this package in publications.
