.onLoad <- function(libname, pkgname) {
  op <- options()
  op.reactome_gsa <- list(
    reactome_gsa.url = "http://193.62.55.4/"
  )
  toset <- !(names(op.reactome_gsa) %in% names(op))
  if(any(toset)) options(op.reactome_gsa[toset])

  invisible()
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("
  +-------------------------------------+
  |                                     |
  |   REACTOME Gene Set Analsis (GSA)   |
  |                                     |
+-+-------------------------------------+-+
|                                         |
|   The ReactomeGSA package provides an   |
|   R interface to the REACTOME Analysis  |
|   Service API.                          |
|                                         |
|   The main feature is to perform gene   |
|   set analysis using the API and view   |
|   the results in REACTOME's web-based   |
|   pathway brower.                       |
|                                         |
+-----------------------------------------+
")
}
