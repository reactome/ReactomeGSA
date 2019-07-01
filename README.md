# ReactomeGSA

The `ReactomeGSA` package is an R client to the `Reactome Analysis System`. This new analysis system supports **multi-species, multi-omics, comparative pathway analyses**.

## Installation

The `ReactomeGSA` package is currently only available on GitHub. Therefore, you need the `devtools` package to install it:

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
