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
    stop("Error: Invalid URL passed. The URL must start with 'http'", call. = FALSE)
  }

  # make sure the url contains a trailing '/'
  if (substr(reactome_url, nchar(reactome_url), nchar(reactome_url)) != "/") {
    reactome_url <- paste0(reactome_url, "/")
  }
  
  # Fixes SEC_E_ILLEGAL_MESSAGE error on Windows 2012
  # This error seems to relate to a problem with Windows' certificate
  # management and the way sonlite::fromJSON acccess the server
  if (substr(reactome_url, 1, 8) == "https://") {
    reactome_url <- tryCatch(
      error = function(cnd) {
        # check if it is really the SEC_E_ILLEGAL_MESSAGE error
        if (grepl("SEC_E_ILLEGAL_MESSAGE", cnd[[1]])) {
          # simply revert to http
          warning("Failed to process certificate. This is generally caused by a problem with Windows and R. Reverting to non-encrypted connection to the ReactomeGSA system.")
          
          # get the non https url
          return(gsub("https://", "http://", reactome_url)[[1]])
        }
        
        # create a nice error if the hostname is unknown
        if (grepl("Could not resolve host", cnd[[1]])) {
          stop("Faild to resolve ", reactome_url, ". Please ensure that the ReactomeGSA URL set under 'options(\"reactome_gsa.url\")' is correct. To change this URL use options(reactome_gsa.url = \"http://gsa.reactome.org\")", call. = FALSE)
        }
        
        # simply fail if there's another issue
        stop(paste0("Failed to connect to ReactomeGSA at '", reactome_url, ": ", cnd[[1]], ". Try reverting the ReactomeGSA url to http://gsa.reactome.org using options(reactome_gsa.url = \"http://gsa.reactome.org\")."), call. = FALSE)
      },
      {
        # a dummy request
        available_methods <- jsonlite::fromJSON(paste0(reactome_url, "0.1/methods"))
        
        # return the default URL
        return(reactome_url)
      }
    )
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