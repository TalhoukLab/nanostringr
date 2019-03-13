nanostringr
===========

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/OVCARE/nanostringr.svg?branch=master)](https://travis-ci.org/OVCARE/nanostringr)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/OVCARE/nanostringr?branch=master&svg=true)](https://ci.appveyor.com/project/OVCARE/nanostringr)
[![Codecov test coverage](https://codecov.io/gh/OVCARE/nanostringr/branch/master/graph/badge.svg)](https://codecov.io/gh/OVCARE/nanostringr?branch=master)
[![CRAN status](https://www.r-pkg.org/badges/version/nanostringr)](https://cran.r-project.org/package=nanostringr)
<!-- badges: end -->

Overview
--------

An R Package for quality assurance checking, normalization and batch effects adjustments of NanoString data suitable for single sample processing. This is the companion R package for the paper ["Single-Patient Molecular Testing with NanoString nCounter Data Using a Reference-Based Strategy for Batch Effect Correction"](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0153844), published in PLOS ONE.


Installation
------------

To install this package, use `devtools`:

``` r
# install.packages(devtools)
devtools::install_github("OVCARE/nanostringr")
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
