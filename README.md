# ReactomeGSA

The `ReactomeGSA` package is an R client to the `Reactome Analysis System`. This new analysis system supports **multi-species, multi-omics, comparative pathway analyses**.

## Getting Help

  * For any questions surrounding the use of ReactomeGSA, please simply post it in our [Q&A Section](https://github.com/reactome/ReactomeGSA/discussions)
  * Should you find a bug in our tool(s) we'd be very grateful if you could report it by posting a [new GitHub Issue](https://github.com/reactome/ReactomeGSA/issues/new)
  * If you encounter any other issues, don't hesitate to contact us at help [at] reactome [dot] org

## Documnetation

The complete usage of the package is described in the [main vignette](https://bioconductor.org/packages/release/bioc/vignettes/ReactomeGSA/inst/doc/using-reactomegsa.html).

If you are interested in processing single-cell RNA-seq data using ReactomeGSA, checkout [this vignette](https://bioconductor.org/packages/release/bioc/vignettes/ReactomeGSA/inst/doc/using-reactomegsa.html)

## Installation

### Bioconductor

The `ReactomeGSA` package is part of Bioconductor since version 3.10. You can find detailed information on the latest stable version on [ReactomeGSA's Bioconductor page](https://doi.org/doi:10.18129/B9.bioc.ReactomeGSA).

To install the latest version from Bioconductor use

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ReactomeGSA")
```

### Latest Version

Bioconductor is updated every 6 months. You can still get the latest version of the ReactomeGSA package directly from GitHub:

```r
# install devtools if needed
if (!require(devtools)) {
  install.packages("devtools")
}

# install the ReactomeGSA package
if (!require(ReactomeGSA)) {
  install_github("reactome/ReactomeGSA")
}
```
