nanostringr
===========

An R Package for quality assurance checking, normalization and batch effects adjustments of NanoString data suitable for single sample processing. This is the companion R package for the paper ["Single-Patient Molecular Testing with NanoString nCounter Data Using a Reference-Based Strategy for Batch Effect Correction"](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0153844), published in PLOS ONE.


Installation
------------

To install this package, use devtools:

``` r
devtools::install_github("OVCARE/nanostringr", build_vignettes = TRUE)
```

Dependencies
------------

Please note that the package [CHL26predictor](https://github.com/tinyheero/CHL26predictor) is needed to run the Hodgkin Lymphoma predictive model.

Overview
--------

To see the full list of exported functions:

``` r
library(nanostringr)
ls("package:nanostringr")
```

A quick overview of the key functions:

-   `NanoStringQC`: Computes quality assurance metrics.
-   `HKnorm`: Performs log (base 2) transformation and normalization to Housekeeping Genes
-   `refMethod`: Performs batch effect correction using the reference-based strategie

A [vignette](http://htmlpreview.github.io/?https://github.com/AlineTalhouk/nanostringr/blob/master/vignettes/my-vignette.html) that reproduces most of the analyses in the paper is included. The vignettes can be accessed in R  using 

``` r
browseVignettes("nanostringr")
```
