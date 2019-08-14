.onLoad <- function(libname, pkgname) {
  op <- options()
  op.reactome_gsa <- list(
    reactome_gsa.url = "http://gsa.reactome.org/"
  )
  toset <- !(names(op.reactome_gsa) %in% names(op))
  if(any(toset)) options(op.reactome_gsa[toset])

  invisible()
}

.onAttach <- function(libname, pkgname) {

}
