#' get_reactome_methods
#'
#' Returns all available analysis methods from the Reactome analysis service.
#'
#' Every method has a type, a scope, and sometimes a list of allowed values. The type (string, int = integer, float) define
#' the expected data type. The \strong{scope} defines at what level the parameter can be set. \emph{dataset} level parameters
#' can be set at the dataset level (using the \code{\link{add_dataset}} function) or at the analysis request level (using
#' \code{\link{set_parameters}}). If these parameters are set at the analysis request level, this overwrites the default
#' value for all datasets. \emph{analysis} and \emph{global} level parameters must only be set at the analysis request level
#' using \code{\link{set_parameters}}. The difference between these two types of parameters is that while \emph{analysis}
#' parameters influence the results, \emph{global} parameters only influence the behaviour of the analysis system (for example
#' whether a Reactome visualization is created).
#'
#' @param print_methods If set to \code{TRUE} (default) a (relatively) nice formatted version of the result is printed.
#' @param print_details If set to \code{TRUE} detailed information about every method, including available parameters and
#'                      description are displayed. This does not affect the data returned if \code{return_result} is \code{TRUE}.
#' @param return_result If set to \code{TRUE}, the result is returned as a data.frame (see below)
#' @param method If set to a method's id, only information for this method will be shown. This is especially useful if
#'               detailed information about a single method should be retrieved. This does not affect the data returned
#'               if \code{return_result} is \code{TRUE}.
#' @param reactome_url URL of the Reactome API Server. Overwrites the URL set in the 'reactome_gsa.url' option.
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
#'
#' # get the details for PADOG
#' get_reactome_methods(print_details = TRUE, method = "PADOG")
get_reactome_methods <- function(print_methods, print_details, return_result, method, reactome_url) {
  # set the default values
  if (missing(print_methods)) print_methods <- TRUE
  if (missing(print_details)) print_details <- FALSE
  if (missing(return_result)) return_result <- FALSE
  if (missing(reactome_url)) reactome_url <- NULL

  reactome_url <- check_reactome_url(reactome_url)

  # fetch the methods
  available_methods <- jsonlite::fromJSON(paste0(reactome_url, "0.1/methods"))

  # print a nice version of the methods
  if (print_methods) {
    for (n_method in seq(to = nrow(available_methods))) {
      # ignore all methods that are not "method" if set
      if (!missing(method) && tolower(available_methods[n_method, "name"]) != tolower(method)) {
        next
      }

      # show a box around the method name if details are shown
      if (print_details) {
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

        cat("\n")
      } else {
        # only print the method name if no details are set
        cat(available_methods[n_method, "name"], ": ", available_methods[n_method, "description"], "\n", sep = "")
      }
    }
  }

  if (return_result) {
    return(available_methods)
  }
}

#' ReactomeGSA supported data types
#'
#' @param print_types If set to \code{TRUE} (default) a (relatively) nice formatted version of the result is printed.
#' @param return_result If set to \code{TRUE}, the result is returned as a data.frame (see below)
#' @param reactome_url URL of the Reactome API Server. Overwrites the URL set in the 'reactome_gsa.url' option.
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
get_reactome_data_types <- function(print_types, return_result, reactome_url) {
  # set the default values
  if (missing(print_types)) print_types <- TRUE
  if (missing(return_result)) return_result <- FALSE
  if (missing(reactome_url)) reactome_url <- NULL

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
