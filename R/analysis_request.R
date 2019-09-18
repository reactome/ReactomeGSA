# ReactomeAnalysisRequest class definition ----

#' ReactomeAnalysisRequest class
#'
#' This class is used to collect all information required to submit an
#' analysis request to the Reactome Analysis System.
#'
#' @slot method character. Name of the method to use
#' @slot request_object list. This slot should not be set manually. It stores the internal
#'                      request representation and should be modified using the classes' functions.
#'                      To add parameters, use \code{\link{set_parameters,ReactomeAnalysisRequest-method}}
#'
#' @return A ReactomeAnalysisRequest object.
#' @import methods
#' @export
#'
#' @examples
#' library(ReactomeGSA.data)
#' library(methods)
#'
#' # create the request method and specify its method
#' request <- ReactomeAnalysisRequest(method = "Camera")
#'
#' # add a dataset to the request
#' data(griss_melanoma_proteomics)
#'
#' request <- add_dataset(request = request,
#'              expression_values = griss_melanoma_proteomics,
#'              name = "Proteomics",
#'              type = "proteomics_int",
#'              comparison_factor = "condition",
#'              comparison_group_1 = "MOCK",
#'              comparison_group_2 = "MCM",
#'              additional_factors = c("cell.type", "patient.id"))
#'
#' # to launch the actual analysis use the perform_reactome_analysis function
ReactomeAnalysisRequest <- setClass("ReactomeAnalysisRequest",
                                       slots = list(method = "character", request_object = "list"),
                                       prototype = prototype(method = "Camera"))

# ---- constructor function

#' @export
ReactomeAnalysisRequest <- function(method) {
  request <- methods::new("ReactomeAnalysisRequest", method = method)
  return(request)
}

#' Check's if a ReactomeAnalysisRequest object is valid
#'
#' @param request The request object to check.
#'
#' @return TRUE if the object is valid or a string with the reason why it is not
checkRequestValidity <- function(request) {
  # the method name must be a string
  if (!methods::is(request@method, "character")) {
    return("'method' must be of class 'character'")
  }

  # method must only contain one element
  if (length(request@method) != 1) {
    return("'method' must be a single string.")
  }

  # make sure that method is not empty
  if (nchar(request@method) < 1) {
    return("'method' must be set")
  }

  if (length(request@request_object) > 0) {
    # TODO: make sure the request object is valid
    if (!all(names(request@request_object) %in% c("parameters", "datasets")) || length(request@request_object) > 2) {
      return("'request_object' must only contain 'parameters' and 'datasets' slots")
    }
  }

  return(TRUE)
}

setValidity("ReactomeAnalysisRequest", checkRequestValidity)

# print and show function ----
#' print - ReactomeAnalysisRequest
#'
#' Shows a \code{\link{ReactomeAnalysisRequest}} object summary.
#'
#' @param object \code{\link{ReactomeAnalysisRequest}}
#' @return The classname of the object
#'
#' @export
#'
#' @examples
#' library(methods)
#'
#' request <- ReactomeAnalysisRequest(method = "Camera")
#' print(request)
#'
#' # add additional parameters
#' request <- set_parameters(request, "max_missing_values" = 0.5)
#' show(request)
setMethod("show", c("object" = "ReactomeAnalysisRequest"), function(object) print(object))

#' print - ReactomeAnalysisRequest
#'
#' Shows a \code{\link{ReactomeAnalysisRequest}} object summary.
#'
#' @param x \code{\link{ReactomeAnalysisRequest}}
#' @param ... Not used
#' @return The classname of the object
#'
#' @export
#'
#' @examples
#' library(methods)
#'
#' request <- ReactomeAnalysisRequest(method = "Camera")
#' print(request)
#'
#' # add additional parameters
#' request <- set_parameters(request, "max_missing_values" = 0.5)
#' show(request)
setMethod("print", c("x" = "ReactomeAnalysisRequest"), function(x, ...) {
  cat("ReactomeAnalysisRequestObject\n")
  cat("  Method = ", x@method, "\n", sep = "")

  if (length(x@request_object) == 0) {
    cat("  No request data stored\n")
  } else {
    # parameters
    if ("parameters" %in% names(x@request_object)) {
      cat("  Parameters:\n")
      for (i in seq_along(rownames(x@request_object[["parameters"]]))) {
        parameter_name <- x@request_object[["parameters"]][i, "name"]
        parameter_value <- x@request_object[["parameters"]][i, "value"]
        cat("  - ", parameter_name, ": ", parameter_value, "\n", sep = "")
      }
    } else {
      cat("  Parameters: none\n")
    }

    # datasets
    if ("datasets" %in% names(x@request_object)) {
      cat("  Datasets:\n")
      for (i in seq_along(x@request_object[["datasets"]][, "name"])) {
        cat("  - ", x@request_object[["datasets"]][i, "name"], " (", x@request_object[["datasets"]][i, "type"], ")\n", sep = "")

        # print the parameters if there are any
        dataset_params <- x@request_object[["datasets"]][i, "parameters"][[1]]

        if (nrow(dataset_params) > 0) {
          for (j in seq(to = nrow(dataset_params))) {
            cat("    ", dataset_params[j, "name"], ": ", dataset_params[j, "value"], "\n", sep = "")
          }
        } else {
          cat("    No parameters set.\n")
        }
      }
    } else {
      cat("  Datasets: none\n")
    }
  }

  return(class(x))
})

# set_parameters function ----

#' set_parameters
#'
#' Sets the analysis parameters for the given \code{\link{ReactomeAnalysisRequest}}. If the
#' parameter is already set, it is overwritten. Use \code{\link{get_reactome_methods}}
#' to get a list of all available parameters for each available method.
#'
#' Both, parameters with the scope "dataset" as well as "analysis" can be set on the analysis
#' level. In this case, these parameters overwrite the system's default values. If a parameter
#' with the scope "dataset" is defined again at the dataset level, this value will overwrite
#' the analysis' scope value for the given dataset.
#'
#' @param request The \code{\link{ReactomeAnalysisRequest}} to set the parameters for.
#' @param ... Any name / value pair to set a parameter (see example). For a complete list of
#'            available parameters use \code{\link{get_reactome_methods}}
#'
#' @return The modified \code{\link{ReactomeAnalysisRequest}} object
#' @export
#'
#' @examples
#' library(methods)
#'
#' # create a request object
#' request <- ReactomeAnalysisRequest(method = "Camera")
#'
#' # add a parameter
#' request <- set_parameters(request, max_missing_values = 0.5, discrete_norm_function = "TMM")
setGeneric("set_parameters", function(request, ...) standardGeneric("set_parameters"))

#' ReactomeAnalysisRequest - set_parameters
#'
#' @inherit set_parameters
setMethod("set_parameters", c("request" = "ReactomeAnalysisRequest"), function(request, ...) {
  # convert into a list
  method_parameters = list(...)

  # create the data.frame to use or use the existing one
  if ("parameters" %in% names(request@request_object)) {
    method_parameters_df <- request@request_object[["parameters"]]
  } else {
    method_parameters_df <- data.frame(stringsAsFactors = F)
  }

  # add all passed parameters
  for (i in seq_along(method_parameters)) {
    parameter_name = names(method_parameters)[i]
    # all values must be strings
    parameter_value = as.character(method_parameters[[parameter_name]])

    # overwrite if it exists
    if ("name" %in% colnames(method_parameters_df) && parameter_name %in% method_parameters_df$name) {
      method_parameters_df$value[method_parameters_df$name == parameter_name] <- parameter_value
    } else {
      # add a new entry
      method_parameters_df <- rbind(method_parameters_df, data.frame("name" = parameter_name, "value" = parameter_value, stringsAsFactors = F))
    }
  }

  # store the new object
  request@request_object[["parameters"]] <- method_parameters_df

  return(request)
})

# add_dataset function ----
#' add_dataset
#'
#' Adds a dataset to the analysis request
#'
#' @param request The request to add the dataset to. Commonly a \code{\link{ReactomeAnalysisRequest}} object.
#' @param expression_values Object containing the expression values of the dataset to add (multiple types supported).
#' @param name character. Name of the dataset. This must be unique within one request.
#' @param type character. The type of the dataset. Get available types using \code{\link{get_reactome_data_types}}
#' @param comparison_factor character. The name of the sample property to use for the main comparison. The sample properties
#'                          are either retrieved from \code{expression_values} or from \code{sample_data}.
#' @param comparison_group_1 character. Name of the first group within \code{comparison_factor} to use for the comparison.
#' @param comparison_group_2 character. Name of the second group within \code{comparison_factor} to use for the comparison.
#' @param sample_data data.frame (optional) data.frame containing the sample metadata of the \code{expression_values}. Depending
#'                    on the object type of \code{expression_values}, this information can also be extracted from there.
#' @param additional_factors vector. Vector of additional sample properties that are used as blocking factors
#'                           (if supported by the chosen analysis method) in the gene set analysis.
#' @param overwrite boolean. If set to \code{TRUE}, datasets with the same name will be overwritten
#' @param ... Additional parameters passed to downstream functions. See the respective documentation of
#'            whether any additional parameters are supported.
#'
#' @return The \code{\link{ReactomeAnalysisRequest}} object with the added dataset
#' @family add_dataset methods
#' @export
#'
#' @examples
#' # create a request using Camera as an analysis
#' library(ReactomeGSA.data)
#' data(griss_melanoma_proteomics)
#' library(methods)
#'
#' my_request <- ReactomeAnalysisRequest(method = "Camera")
#'
#' # since the expression_values object is a limma EList object, the sample_data is
#' # retrieved from there
#'
#' # add the dataset
#' my_request <- add_dataset(request = my_request,
#'                           expression_values = griss_melanoma_proteomics,
#'                           name = "Proteomics",
#'                           type = "proteomics_int",
#'                           comparison_factor = "condition",
#'                           comparison_group_1 = "MOCK",
#'                           comparison_group_2 = "MCM",
#'                           additional_factors = c("cell.type", "patient.id"))
setGeneric("add_dataset",
           function(request, expression_values, name, type, comparison_factor,
                    comparison_group_1, comparison_group_2, sample_data = NULL, additional_factors, overwrite, ...) {
             standardGeneric("add_dataset")
            })

#' add_dataset - EList
#'
#' @param request ReactomeAnalysisRequest.
#' @param expression_values EList. Here, the sample_data is automaticall extracted from the \code{expression_values} object unless \code{sample_data}
#'                          is specified as well.
#'
#' @family add_dataset methods
#' @inherit add_dataset
setMethod("add_dataset", c("request" = "ReactomeAnalysisRequest", "expression_values" = "EList"),
          function(request, expression_values, name, type, comparison_factor,
                   comparison_group_1, comparison_group_2, sample_data = NULL, additional_factors = NULL, overwrite = FALSE, ...) {
            # make sure the EList contains sample data
            if (is.null(sample_data) && !"samples" %in% names(expression_values)) {
              stop("Error: Missing 'samples' in expression_values EList object")
            }

            # sample.data overwrites the stored one
            if (is.null(sample_data)) {
              sample_data <- expression_values$samples
            }

            expression_values <- data.frame(expression_values$E)

            # call the data.frame implementation
            return(add_dataset(request = request, name = name, expression_values = expression_values, sample_data = sample_data,
                               type = type, comparison_factor = comparison_factor, comparison_group_1 = comparison_group_1,
                               comparison_group_2 = comparison_group_2, additional_factors = additional_factors, overwrite = overwrite,
                               ...))
          })

#' add_dataset - DGEList
#'
#' @param request ReactomeAnalysisRequest.
#' @param expression_values DGEList Here, the sample_data is automaticall extracted from the \code{expression_values} object unless \code{sample_data}
#'                          is specified as well.
#'
#' @family add_dataset methods
#' @inherit add_dataset
setMethod("add_dataset", c("request" = "ReactomeAnalysisRequest", "expression_values" = "DGEList"),
          function(request, expression_values, name, type, comparison_factor,
                   comparison_group_1, comparison_group_2, sample_data = NULL, additional_factors = NULL,
                   overwrite = FALSE, ...) {
            # make sure the EList contains sample data
            if (is.null(sample_data) && !"samples" %in% names(expression_values)) {
              stop("Error: Missing 'samples' in expression_values EList object")
            }

            # sample.data overwrites the stored one
            if (is.null(sample_data)) {
              sample_data <- expression_values$samples
            }

            expression_values <- data.frame(expression_values$counts)

            # call the data.frame implementation
            return(add_dataset(request = request, name = name, expression_values = expression_values, sample_data = sample_data,
                               type = type, comparison_factor = comparison_factor, comparison_group_1 = comparison_group_1,
                               comparison_group_2 = comparison_group_2, additional_factors = additional_factors,
                               overwrite = overwrite, ...))
          })

#' add_dataset - ExpressionSet
#'
#' @param request ReactomeAnalysisRequest.
#' @param expression_values ExpressionSet. Here, the sample_data is automaticall extracted from the \code{expression_values} object unless \code{sample_data}
#'                          is specified as well.
#'
#' @family add_dataset methods
#' @inherit add_dataset
setMethod("add_dataset", c("request" = "ReactomeAnalysisRequest", "expression_values" = "ExpressionSet"),
          function(request, expression_values, name, type, comparison_factor,
                   comparison_group_1, comparison_group_2, sample_data = NULL, additional_factors = NULL,
                   overwrite = FALSE, ...) {
            # make sure the ExpressionSet contains sample data
            has_pdata <- !is.null(Biobase::pData(expression_values)) && nrow(Biobase::pData(expression_values)) > 0
            if (is.null(sample_data) && !has_pdata) {
              stop("Error: Missing 'samples' in expression_values ExpressionSet object")
            }

            # sample.data overwrites the stored one
            if (is.null(sample_data)) {
              sample_data <- Biobase::pData(expression_values)
            }

            expression_values <- data.frame( Biobase::exprs(expression_values) )

            # call the data.frame implementation
            return(add_dataset(request = request, name = name, expression_values = expression_values, sample_data = sample_data,
                               type = type, comparison_factor = comparison_factor, comparison_group_1 = comparison_group_1,
                               comparison_group_2 = comparison_group_2, additional_factors = additional_factors,
                               overwrite = overwrite, ...))

          })

#' add_dataset - data.frame
#'
#' @param request ReactomeAnalysisRequest.
#' @param expression_values data.frame. In this case, the \code{sample_data} must be set.
#'
#' @family add_dataset methods
#' @inherit add_dataset
setMethod("add_dataset", c("request" = "ReactomeAnalysisRequest", "expression_values" = "data.frame"),
          function(request, expression_values, name, type, comparison_factor,
                   comparison_group_1, comparison_group_2, sample_data, additional_factors, overwrite,
                   ...) {
            # set the default value for "overwrite"
            if (missing(overwrite)) overwrite <- FALSE

            # get the dataset_df object or create a new one if none exists
            if ("datasets" %in% names(request@request_object)) {
              dataset_df <- request@request_object[["datasets"]]
            } else {
              dataset_df <- data.frame(stringsAsFactors = F)
            }

            # check if the dataset already exists
            if (!overwrite && "name" %in% colnames(dataset_df) && name %in% dataset_df$name) {
              stop("Error: The request already contains a dataset named '", name, "'")
            }

            # make sure the sample_data contains information
            if (is.null(sample_data) || nrow(sample_data) < 1) {
              stop("Error: Missing required sample_data or no entries (= rows) found in sample_data data.frame")
            }

            # make sure the sample_data fits the expression_values
            if (nrow(sample_data) != ncol(expression_values)) {
              stop("Error: sample_data has ", nrow(sample_data), " rows while expression_data contains ", ncol(expression_values), " columns. Both must match")
            }

            # make sure the comparison_factor is present in sample_data
            if (!comparison_factor %in% colnames(sample_data)) {
              stop("Error: comparison_factor '", comparison_factor, "' is not defined (not in colnames) in the passed sample_data")
            }

            # use the sample_data rownames as sample names
            sample_names <- rownames(sample_data)
            comparison_group <- as.character(sample_data[, comparison_factor])

            # make sure the comparison groups are annotated
            if (!comparison_group_1 %in% comparison_group) {
              stop("Error: No comparison_group_1 samples found. No samples annotated as '", comparison_group_1, "' in '", comparison_factor, "'")
            }
            if (!comparison_group_2 %in% comparison_group) {
              stop("Error: No comparison_group_2 samples found. No samples annotated as '", comparison_group_2, "' in '", comparison_factor, "'")
            }

            # delete the dataset if it exists
            if (name %in% dataset_df$name) {
              message("Overwriting ", name, "...")
              dataset_df <- dataset_df[dataset_df$name != name, ]
            }

            # create the list representing the design
            design <- list(
              "analysisGroup" = comparison_group,
              "comparison" = list(
                "group1" = comparison_group_1,
                "group2" = comparison_group_2
              ),
              "samples" = sample_names)

            # add additional factors if set
            if (!missing(additional_factors)) {
              if (methods::is(additional_factors, "list") || methods::is(additional_factors, "character")) {
                for (factor_name in additional_factors) {
                  if (!factor_name %in% colnames(sample_data)) {
                    stop("Error: Failed to find additional_factor '", factor_name, "' in the passed sample_data (not in colnames)")
                  }
                  design[[factor_name]] <- as.character(sample_data[, factor_name])
                }
              }
            }

            # add all additional parameters (...) as dataset-level parameters
            dataset_parameters <- data.frame(stringsAsFactors = F)
            user_params <- list(...)
            for (param_name in names(user_params)) {
              dataset_parameters <- rbind(dataset_parameters, data.frame(name = param_name,
                                                                         value = as.character(user_params[[param_name]]),
                                                                         stringsAsFactors = F))
            }

            # Create the dataset object as a data.frame row
            dataset_row <- data.frame(
              "data" = data_frame_as_string(expression_values),
              "name" = name,
              "type" = as.character(type),
              stringsAsFactors = F
            )

            dataset_row$design <- list(design)
            dataset_row$parameters <- list(dataset_parameters)

            # save the new dataset
            request@request_object[["datasets"]] <- rbind(dataset_df, dataset_row)

            return(request)
          })

# remove_dataset ----
#' remove_dataset
#'
#' Remove the dataset from the \code{\link{ReactomeAnalysisRequest}} object.
#'
#' @param x The \code{\link{ReactomeAnalysisRequest}} to remove the dataset from
#' @param dataset_name character The dataset's name
#'
#' @return The updated \code{\link{ReactomeAnalysisRequest}}
#' @export
#'
#' @examples
#' # create a request using Camera as an analysis
#' library(ReactomeGSA.data)
#' data(griss_melanoma_proteomics)
#' library(methods)
#'
#' my_request <- ReactomeAnalysisRequest(method = "Camera")
#'
#' # since the expression_values object is a limma EList object, the sample_data is
#' # retrieved from there
#'
#' # add the dataset
#' my_request <- add_dataset(request = my_request,
#'                           expression_values = griss_melanoma_proteomics,
#'                           name = "Proteomics",
#'                           type = "proteomics_int",
#'                           comparison_factor = "condition",
#'                           comparison_group_1 = "MOCK",
#'                           comparison_group_2 = "MCM",
#'                           additional_factors = c("cell.type", "patient.id"))
#'
#' # remove the dataset again
#' my_request <- remove_dataset(x = my_request, dataset_name = "Proteomics)
setGeneric("remove_dataset", function(x, dataset_name) standardGeneric("remove_dataset"))

#' remove_dataset - ReactomeAnalysisRequest
#'
#' @inherit remove_dataset
setMethod("remove_dataset", c("x" = "ReactomeAnalysisRequest"), function(x, dataset_name) {
  # make sure the dataset exists
  if (!"datasets" %in% names(x@request_object)) {
    stop("Error: ReactomeAnalysisRequest object does not contain any datasets.")
  }
  if (!dataset_name %in% x@request_object[["datasets"]][, "name"]) {
    stop("Error: Dataset '", dataset_name, "' does not exist in the passed ReactomeAnalysisRequest object.")
  }

  # get the dataset's index
  dataset_index <- which(x@request_object[["datasets"]][, "name"] == dataset_name)

  # remove the dataset from the table
  x@request_object[["datasets"]] <- x@request_object[["datasets"]][dataset_index * -1, ]

  return(x)
})

# to JSON ----
setGeneric("toJSON", function(x, ...) standardGeneric("toJSON"))

setMethod("toJSON", c("x" = "ReactomeAnalysisRequest"), function(x, ...) {
  # make sure the object is valid
  validity_check <- methods::validObject(x)
  if (validity_check != TRUE) {
    stop("Error: Invalid ReactomeAnalysisRequest object passed. ", validity_check)
  }

  # create the JSON object
  json_request <- jsonlite::toJSON(x@request_object, pretty=FALSE)

  # add the request method using string replacement
  json_request <- paste0("{\"methodName\": \"",  x@method, "\", ", substr(json_request, 2, nchar(json_request)))

  # fix the comparison groups (stored as list, must be single values)
  json_request <- gsub("group([12])\":\\[\"([^\"]*)\"\\]", "group\\1\":\"\\2\"", json_request)

  return(json_request)
})

# set_method -----

#' set_method
#'
#' Set the analysis method used by the \code{\link{ReactomeAnalysisRequest}}
#'
#' @param request The \code{\link{ReactomeAnalysisRequest}} to adjust
#' @param method The name of the method to use. Use \code{\link{get_reactome_methods}} to
#'               retrieve all available methods
#' @param ... Additional parameters passed to specific implementations
#'
#' @return The \code{\link{ReactomeAnalysisRequest}} with the adapted method
#' @export
#'
#' @examples
#' # create a request using Camera as an analysis
#' data(griss_melanoma_proteomics)
#' library(methods)
#'
#' my_request <- ReactomeAnalysisRequest(method = "Camera")
#'
#' print(my_request)
#'
#' # change the method to ssGSEA
#' my_request <- set_method(my_request, "ssGSEA")
#'
#' print(my_request)
setGeneric("set_method", function(request, method, ...) standardGeneric("set_method"))

#' set_method - ReactomeAnalysisRequest
#'
#' @inherit set_method
setMethod("set_method", c("request" = "ReactomeAnalysisRequest"), function(request, method, ...) {
  # only adapt the method
  request@method <- method

  return(request)
})
