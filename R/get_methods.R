#' Returns all available analysis methods from the REACTOME service
#'
#' @param print_methods If set to \code{TRUE} (default) a (relatively) nice formatted version of the result is printed.
#' @param return_result If set to \code{TRUE}, the result is returned as a data.frame (see below)
#' @param reactome_url URL of the REACTOME API Server. Overwrites the URL set in the 'reactome_gsa.url' option.
#'                     Specific ports can be set using the standard URL specification (for example http://your.service:1234)
#' @return If \code{return_result} is set to \code{TRUE}, a data.frame with one row per method. Each method has a name, description, and
#'         (optional) a list of parameters. Parameters again have a name, type, and description.
#' @author Johannes Griss
#' @export
#' @importFrom methods is
#'
#' @family Reactome Service functions
#'
#' @examples
#' # retrieve the available methods only in an object
#' available_methods <- get_reactome_methods(print_methods = FALSE, return_result = TRUE)
#'
#' # print all method names
#' available_methods$name
#'
#' # list all parameters for the first method
#' first_method_parameters <- available_methods[1, "parameters"]
#' first_method_parameters
#'
#' # simply print the available methods
#' get_reactome_methods()
get_reactome_methods <- function(print_methods=TRUE, return_result=FALSE, reactome_url=NULL) {
  reactome_url <- check_reactome_url(reactome_url)

  # fetch the methods
  available_methods <- jsonlite::fromJSON(paste0(reactome_url, "0.1/methods"))

  # print a nice version of the methods
  if (print_methods) {
    for (n_method in seq(to = nrow(available_methods))) {
      method_name_length = nchar(available_methods[n_method, "name"])
      cat("+-", rep("-", method_name_length), "--+\n", sep = "")
      cat("| ", available_methods[n_method, "name"], ": |\n", sep = "")
      cat("+-", rep("-", method_name_length), "--+\n\n", sep = "")
      cat(available_methods[n_method, "description"], "\n\n")
      cat("Parameters:\n")
      if ("parameters" %in% colnames(available_methods) && is(available_methods[n_method, "parameters"], "list")) {
        parameter_df = available_methods[n_method, "parameters"][[1]]
        for (n_param in seq(to = nrow(parameter_df))) {
          cat("- ", parameter_df[n_param, "name"], "\n", sep = "")
          cat("  Type: ", parameter_df[n_param, "type"], "\n", sep = "")
          cat("  Default: ", parameter_df[n_param, "default"], "\n", sep = "")
          cat("  Scope: ", parameter_df[n_param, "scope"], "\n", sep = "")
          if (!is.null(parameter_df[n_param, "values"]) && parameter_df[n_param, "values"] != "NULL") {
            cat("  Allowed values: ", paste0(parameter_df[n_param, "values"][[1]], collapse = ", "), "\n", sep = "")
          }
          cat("  Description: ", parameter_df[n_param, "description"], "\n", sep = "")
        }
      }
    }
  }

  if (return_result) {
    return(available_methods)
  }
}

#' REACTOMEgsa supported data types
#'
#' @param print_types If set to \code{TRUE} (default) a (relatively) nice formatted version of the result is printed.
#' @param return_result If set to \code{TRUE}, the result is returned as a data.frame (see below)
#' @param reactome_url URL of the REACTOME API Server. Overwrites the URL set in the 'reactome_gsa.url' option.
#'                     Specific ports can be set using the standard URL specification (for example http://your.service:1234)
#' @return A \code{data.frame} containing one row per data type with its \code{id} and \code{description}.
#' @author Johannes Griss
#' @export
#' @family Reactome Service functions
#'
#' @examples
#' # retrieve the avialable data types
#' available_types <- get_reactome_data_types(print_types = FALSE,  return_result = TRUE)
#'
#' # print all data type ids
#' available_types$id
#'
#' # simply print the available methods
#' get_reactome_data_types()
get_reactome_data_types <- function(print_types=TRUE, return_result=FALSE, reactome_url=NULL) {
  reactome_url <- check_reactome_url(reactome_url)

  # fetch the types
  available_types <- jsonlite::fromJSON(paste0(reactome_url, "0.1/types"))

  # print a nice version of the methods
  if (print_types) {
    for (n_type in seq(to = nrow(available_types))) {
      cat(available_types[n_type, "id"], ":\n", sep = "")
      cat("  ", available_types[n_type, "name"], "\n  ", available_types[n_type, "description"], "\n")
    }
  }

  if (return_result) {
    return(available_types)
  }
}
