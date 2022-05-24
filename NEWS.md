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
