# ReactomeGSA

The `ReactomeGSA` package is an R client to the `Reactome Analysis System`. This new analysis system supports **multi-species, multi-omics, comparative pathway analyses**.

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

## Help

The complete usage of the package is described in the [main vignette](./vignettes/using-reactomegsa.Rmd).
