# nanostringr 0.6.1

* allow normalization of pools from same CodeSet as the reference pools, e.g. when we want to normalize cross-site replicates. Add parameter `same_codeset` to `normalize_pools()` so we can use a different file name pattern to match for these pools.

# nanostringr 0.6.0

* Replace magrittr pipe with native pipe, subsequently increase minimum R version to 4.2
* Update `read_rcc()` to replace gene names PD1 with PD-1 and PDL1 with PD-L1

# nanostringr 0.5.0

* Fix `normalize_pools()` and `normalize_random()` functions so they actually work
* Remove `magrittr` from Suggests
* Remove deprecated usage of `.data`

# nanostringr 0.4.2

* Fix bug in `NanoStringQC()` where the signal to noise ratio input threshold was not being compared against the observed SNR because the parameter name was the same as the column name. Hence the comparison was against itself, and `normFlag` only checked the percent of genes detected.

# nanostringr 0.4.1

* Fix warning regarding old-style CITATION
* Add ORCID to package authors
* Use GitHub Actions for pkgdown site
* use stricter R CMD check standard workflow
* decrease package dependencies: removed `epiR`, `forcats`, `tibble` from Imports, moved `magrittr` to Suggests

# nanostringr 0.4.0

* new `normalize_pools()` function to use common pool samples to correct for batch effects
when normalizing
* new `normalize_random()` function to use randomly selected samples as the reference for `refMethod()`

# nanostringr 0.3.0

* fix calculation of genes detected so that the limit of detection is compared in parallel to each sample's counts instead of being recycled, #20

* use GitHub Actions for CI, replacing Travis and Appveyor

* add gene label "PC_1" for checking smallest positive control

* more robust extraction of numeric concentrations from PC gene labels

# nanostringr 0.2.0

* use RCC file names in parsed data of `read_rcc()`. Also rename gene name CD3E to CD3e for compatibility purposes

* update roxygen

* remove Rplots.pdf generated from tests, removed deprecated `context()`

# nanostringr 0.1.4

* reduce package dependencies from Imports

* remove old packages from Suggests

* performance gains in `refMethod()`

* `raw` parameter in `NanoStringQC()` now accepts tibbles

# nanostringr 0.1.3

* check for presence of "Positive" and "Negative" genes in `NanoStringQC()`

* remove most code lints

* unit tests for `read_rcc()` and internal functions `check_colnames()`, `check_genes()`

* fix bug using `NanoStringQC()` for single sample data

# nanostringr 0.1.2

* sort sample names in HLD and OVD cohorts by numeric order, not lexicographic order

* tidy up code in `HKnorm()`

* use tidy evaluation in `NanostringQC()`

* use `vapply()` over `sapply()` for input checking functions

# nanostringr 0.1.1

* update citation to link to PLOS paper

# nanostringr 0.1.0

* New submission to CRAN accepted on March 15, 2019
