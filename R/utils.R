#' check_reactome_url
#'
#' Makes sure the passed URL is valid. If not URL is passed, the one
#' stored in the options is retrieved
#'
#' @param reactome_url character The URL to test. If \code{NULL} the URL
#'                     is retrieved from the options.
#' @return character The potentially cleaned / retrieved URL with a trailing "/"
check_reactome_url <- function(reactome_url) {
  if (is.null(reactome_url)) {
    reactome_url <- getOption("reactome_gsa.url")
  }

  if (substr(reactome_url, 1, 4) != "http") {
    stop("Error: Invalid URL passed. The URL must start with 'http'")
  }

  # make sure the url contains a trailing '/'
  if (substr(reactome_url, nchar(reactome_url), nchar(reactome_url)) != "/") {
    reactome_url <- paste0(reactome_url, "/")
  }

  return(reactome_url)
}

#' Converts a data.frame to a string representation
#'
#' A data.frame is converted into a single string using
#' `\\t` (the characters, not tab) as field delimiter and
#' `\\n` (the characters, not newline) as line delimiter
#'
#' @param data The data.frame to convert
#'
#' @return A string representing the passed data.frame
data_frame_as_string <- function(data) {
  if (!methods::is(data, "data.frame")) {
    stop("Error: Only data.frame objects can be converted to string representation")
  }

  if (nrow(data) > 10000) {
    message("Converting expression data to string... (This may take a moment)")
  }

  text.connection <- textConnection(NULL, "w")
  utils::write.table(data, text.connection, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
  text.result <- textConnectionValue(text.connection)
  close(text.connection)

  text.result <- paste0(text.result, collapse = "\n")

  if (nrow(data) > 10000) {
    message("Conversion complete")
  }

  return(text.result)
}
