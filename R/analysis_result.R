# ReactomeAnalysisResult class definition --------------------------------------------

#' ReactomeAnalysisResult class
#'
#' A ReactomeAnalysisResult object contains the pathway analysis results of all
#' submitted datasets at once.
#'
#' This class represents a result retrieved from the Reactome Analysis Service. It is returned
#' by \code{\link{get_reactome_analysis_result}} and its wrapper \code{\link{perform_reactome_analysis}}.
#' Generally, object of this class should not be created manually.
#'
#' @section Methods:
#'
#' \bold{\code{\link{names}}}:
#' Retrieves the names of all datasets in the result object
#'
#' \bold{\code{\link{result_types}}}:
#' Retrieves the available result types
#'
#' \bold{\code{\link{pathways}}}:
#' Merges the pathway results of all analysed datasets.
#'
#' \bold{\code{\link{get_result}}}:
#' Retrieve a specific result as data.frame
#'
#' \bold{\code{\link{reactome_links}}}:
#' Displays / retrieves the URLs to the available visualizations in Reactome's pathway browser.
#'
#' \bold{\code{\link{open_reactome}}}:
#' Opens the specified Reactome visualization in the system's default browser.
#'
#' @slot reactome_release The Reactome version used to create this result.
#' @slot mappings Stores the mapping results that were generated for this analysis.
#' @slot results A named list containing the actual analysis results for every dataset
#'               and possibly combined results as well.
#' @slot reactome_links Links pointing to reactome results as a list.
#'
#' @return A ReactomeAnalysisResult object.
#' @export
#' @import methods
#'
#' @examples
#' # load an example result object
#' library(ReactomeGSA.data)
#' data(griss_melanoma_result)
#'
#' # retrieve the names of all datasets in the result
#' names(griss_melanoma_result)
#'
#' # get the combined pathway result
#' pathway_result <- pathways(griss_melanoma_result)
#'
#' # check which result types are available
#' result_types(griss_melanoma_result)
#'
#' # get the fold changes for the first dataset
#' first_dataset_name <- names(griss_melanoma_result)[1]
#'
#' first_fc <- get_result(griss_melanoma_result, "fold_changes", first_dataset_name)
ReactomeAnalysisResult <- setClass("ReactomeAnalysisResult",
                                    slots = list(reactome_release = "character", mappings = "list", results = "list", reactome_links = "list"))

# create request function ---------------------------------------------

#' Convert the Reactome JSON result to a ReactomeAnalysisResult object
#'
#' @param reactome_result The JSON result already converted to R objects (name list)
#' @import methods
#'
#' @return A \code{\link{ReactomeAnalysisResult}} object
convert_reactome_result <- function(reactome_result) {
  loadNamespace("methods")

  # make sure the object has the required fields
  required_fields <- c("release", "mappings", "results")
  for (field in required_fields) {
    if (!field %in% names(reactome_result)) {
      stop("Error: Invalid Reactome result JSON object found. Missing required field '", field, "'")
    }
  }

  # convert the mapping result
  mapping_result <- list()

  if ("mappings" %in% names(reactome_result)) {
    for (i in seq(to = nrow(reactome_result[["mappings"]]))) {
      mapping_result[[reactome_result[["mappings"]][i, "identifier"]]] <- reactome_result[["mappings"]][i, "mapped_to"]
    }
  }

  # convert the results
  results <- list()

  for (i in seq(to = nrow(reactome_result[["results"]]))) {
    cur_result <- reactome_result[["results"]][i, ]

    cur_name <- cur_result[1, "name"]
    pathway_string <- cur_result[1, "pathways"]
    fc_string <- cur_result[1, "fold_changes"]

    # represent one result as a list object
    results[[cur_name]] <- list()

    if (!is.null(pathway_string) && nchar(pathway_string) > 0) {
      results[[cur_name]][["pathways"]] <- utils::read.csv(text = pathway_string, sep = "\t")
    }

    if (!is.na(fc_string) && nchar(fc_string) > 0) {
      results[[cur_name]][["fold_changes"]] <- utils::read.csv(text = fc_string, sep = "\t")
    }
  }

  # add the reactome links
  reactome_links = list()
  if ("reactome_links" %in% names(reactome_result) && length(reactome_result[["reactome_links"]]) > 0) {
    for (i in seq(to = nrow(reactome_result[["reactome_links"]]))) {
      reactome_links[[i]] <- c(url = reactome_result[["reactome_links"]][i, "url"],
                               name = reactome_result[["reactome_links"]][i, "name"],
                               description = reactome_result[["reactome_links"]][i, "description"])
    }
  }

  # create the result object
  result_object <- methods::new("ReactomeAnalysisResult",
                       reactome_release = reactome_result[["release"]],
                       mappings = mapping_result,
                       results = results,
                       reactome_links = reactome_links)

  return(result_object)
}

# print and show function ----
#' show - ReactomeAnalysisResult
#'
#' Displays basic information about the \code{\link{ReactomeAnalysisResult}}
#' object.
#'
#' @param object ReactomeAnalysisResult.
#' @return character classname of the object
#' @export
#'
#' @examples
#' library(ReactomeGSA.data)
#' data(griss_melanoma_result)
#'
#' show(griss_melanoma_result)
setMethod("show", c("object" = "ReactomeAnalysisResult"), function(object) print(object))

#' print - ReactomeAnalysisResult
#'
#' Displays basic information about the \code{\link{ReactomeAnalysisResult}}
#' object.
#'
#' @param x ReactomeAnalysisResult.
#' @param ... Not used
#' @return character classname of the object
#'
#' @export
#'
#' @examples
#' library(ReactomeGSA.data)
#' data(griss_melanoma_result)
#'
#' print(griss_melanoma_result)
setMethod("print", c("x" = "ReactomeAnalysisResult"), function(x, ...) {
  cat("ReactomeAnalysisResult object\n")
  cat("  Reactome Release: ", x@reactome_release, "\n", sep = "")

  if (length(x@results) < 0) {
    cat("  - no results -\n")
  } else {
    cat("  Results:\n")
  }

  for (result_name in names(x@results)) {
    cat("  - ", result_name, ":\n", sep = "")

    if ("pathways" %in% names(x@results[[result_name]])) {
      cat("    ", nrow(x@results[[result_name]][["pathways"]]), " pathways\n", sep = "")
    }

    if ("fold_changes" %in% names(x@results[[result_name]])) {
      cat("    ", nrow(x@results[[result_name]][["fold_changes"]]), " fold changes for genes\n", sep = "")
    }
  }

  # Display the Reactome visualizations
  if (length(x@reactome_links) == 0) {
    cat("  No Reactome visualizations available\n")
  } else {
    cat("  Reactome visualizations:\n")
    for (link in x@reactome_links) {
      cat("  - ", link["name"], "\n", sep = "")
    }
  }

  return(class(x))
})

# names ----------------------------------------------------------------

#' ReactomeAnalysisResult - names
#'
#' Retrieves the names of the contained datasets within an \code{\link{ReactomeAnalysisResult}}
#' object.
#'
#' @param x ReactomeAnalysisResult.
#'
#' @return character vector with the names of the contained datasets
#' @export
#' @family ReactomeAnalysisResult functions
#'
#' @examples
#' # load an example result object
#' library(ReactomeGSA.data)
#' data(griss_melanoma_result)
#'
#' # get the names of the available datasets
#' names(griss_melanoma_result)
setMethod("names", c("x" = "ReactomeAnalysisResult"), function(x) {
  return(names(x@results))
})

# result_types --------------------------------------------------------
#' result_types
#'
#' Retrieves the available result types for the \code{\link{ReactomeAnalysisResult}} object. Currently,
#' the Reactome Analysis System supports \code{pathways} and gene level \code{fold_changes}
#' as result types. Not all analysis methods return both data types though.
#' Use the \code{names} function to find out which datasets are available
#' in the result object.
#'
#' @param x ReactomeAnalysisResult.
#'
#' @return A cacharacter vector of result types.
#' @export
#' @family ReactomeAnalysisResult functions
#'
#' @examples
#' # load an example result object
#' library(ReactomeGSA.data)
#' data(griss_melanoma_result)
#'
#' # get the available result types
#' result_types(griss_melanoma_result)
setGeneric("result_types", function(x) standardGeneric("result_types"))

#' ReactomeAnalysisResult - result_types
#'
#' @inherit result_types
setMethod("result_types", c("x" = "ReactomeAnalysisResult"), function(x) {
  result_types <- c()

  for (result_name in names(x@results)) {
    result_types <- c(result_types, names(x@results[[result_name]]))
  }

  return(unique(result_types))
})

# get_result ------------------------------------------------------

#' get_result
#'
#' Retrieves a result from a \code{\link{ReactomeAnalysisResult}} object.
#'
#' @param x ReactomeAnalysisResult.
#' @param type the type of result. Use \code{\link{result_types}}
#'             to retrieve all available types.
#' @param name the name of the result. Use \code{\link{names}}
#'             to retrieve all available results.
#'
#' @return A \code{data.frame} containing the respective result.
#' @export
#' @family ReactomeAnalysisResult functions
#'
#' @examples
#' # load an example result object
#' library(ReactomeGSA.data)
#' data(griss_melanoma_result)
#'
#' # get the available result types
#' result_types(griss_melanoma_result)
#'
#' # get the dataset names
#' names(griss_melanoma_result)
#'
#' # get the fold_changes for the first dataset
#' prot_fc <- get_result(griss_melanoma_result, type = "fold_changes", name = "proteomics")
#'
#' head(prot_fc)
setGeneric("get_result", function(x, type, name) standardGeneric("get_result"))

#' ReactomeAnalysisResult - get_result
#'
#' @inherit get_result
setMethod("get_result", c("x" = "ReactomeAnalysisResult"), function(x, type, name) {
  if (!name %in% names(x@results)) {
    stop("Error: Failed to find '", name, "' in ReactomeAnalysisResult object")
  }
  if (!type %in% names(x@results[[name]])) {
    stop("Error: Result '", name, "' does not have a '", type, "' result")
  }

  return(x@results[[name]][[type]])
})

# pathways -----------------------------------------------------------------------------
#' pathways
#'
#' Combines and returns the pathways of all analysed datasets.
#'
#' @param x ReactomeAnalysisResult.
#' @param ... Additional parameters for specific implementations.
#'
#' @return A \code{data.frame} containing all merged pathways.
#' @export
#' @family ReactomeAnalysisResult functions
#'
#' @examples
#' # load an example result
#' library(ReactomeGSA.data)
#' data(griss_melanoma_result)
#'
#' # get the combined pathway result
#' pathway_result <- pathways(griss_melanoma_result)
#'
#' head(pathway_result)
setGeneric("pathways", function(x, ...) standardGeneric("pathways"))

#' ReactomeAnalysisResult - pathways
#'
#' @param p Minimum p-value to accept a pathway as significantly regulated. Default is 0.01.
#' @param order_by Name of the dataset to sort the result list by. By default, the
#'                 results are sorted based on the first dataset.
#'
#' @inherit pathways
setMethod("pathways", c("x" = "ReactomeAnalysisResult"), function(x, p, order_by, ...) {
  combined_result <- data.frame()

  # set the default value for p
  if (missing(p)) p <- 0.01

  params <- list(...)

  # merge the pathways
  for (result_name in names(x@results)) {
    if ("pathways" %in% names(x@results[[result_name]])) {
      tmp_data <- x@results[[result_name]][["pathways"]]

      # use the pathway as rowname
      rownames(tmp_data) <- tmp_data[, "Pathway"]
      tmp_data[, "Pathway"] <- NULL

      # add the "sig" column for GSA based approaches
      if ("FDR" %in% colnames(tmp_data)) {
        tmp_data$sig <- tmp_data$FDR <= p
      }

      # add a suffix to the data
      colnames(tmp_data) <- paste0(colnames(tmp_data), ".", result_name)

      combined_result <- merge(combined_result, tmp_data, by.x = 0, by.y = 0,
                               all.x = T, all.y = T, suffixes = c("", paste0(".", result_name)))

      # remove the redundant name column - if present
      new_name_column <- paste0("Name.", result_name)

      if (new_name_column %in% colnames(combined_result)) {
        if ("Name" %in% colnames(combined_result)) {
          missing_names <- is.na(combined_result$Name)
          combined_result$Name[missing_names] <- as.character(combined_result[missing_names, new_name_column])
        } else {
          combined_result$Name <- as.character( combined_result[, new_name_column] )
        }
      }

      combined_result[, new_name_column] <- NULL

      # use the "Row.names" as row.names
      rownames(combined_result) <- combined_result[, "Row.names"]
      combined_result[, "Row.names"] <- NULL
    }
  }

  # sort the result if the parameter is set
  if (!missing(order_by) && !order_by %in% names(x@results)) {
    warning("Warning: order_by dataset '", order_by, "' does not exist. Ignoring parameter.")
    order_by <- names(x@results)[1]
  }

  # if not set, sort according to the first dataset
  if (missing(order_by) || is.null(order_by)) {
    order_by <- names(x@results)[1]
  }

  # use the PValue of that dataset
  order_by <- paste0("PValue.", order_by)[[1]]

  # make sure this column is present - which is not the case in GSVA results
  if (order_by %in% colnames(combined_result)) {
    combined_result <- combined_result[order(combined_result[, order_by]), ]
  }

  # move the "name" as the second column
  if ("Name" %in% colnames(combined_result)) {
    name_col <- which(colnames(combined_result) == "Name")
    non_name_cols <- which(colnames(combined_result) != "Name")
    combined_result <- combined_result[, c(name_col, non_name_cols)]
  }

  return(combined_result)
})

# visualizations ------------------------------------------
#' reactome_links
#'
#' Displays detailed information about the result visualizations in Reactome.
#'
#' @param x ReactomeAnalysisResult.
#' @param ... Additional parameters for specific implementations.
#' @return If \code{return_result} is set to \code{TRUE}, a vector of the available visualizations.
#' @export
#' @family ReactomeAnalysisResult functions
#'
#' @examples
#' # Note: This function only works with a newly created result
#' # since the visualization links only stay active for 7 days
#'
#' # load an example result
#' library(ReactomeGSA.data)
#' data(griss_melanoma_result)
#'
#' # get the reactome link - this does only work
#' # with new results
#' reactome_links(griss_melanoma_result)
setGeneric("reactome_links", function(x, ...) standardGeneric("reactome_links"))

#' ReactomeAnalysisResult - reactome_links
#'
#' @param print_result If set to \code{FALSE} the links are not printed to the console.
#' @param return_result If \code{TRUE} the available visualizations are returned as a list containing
#'                      named vectors for every visualization. These vectors' have a \code{url}, \code{name},
#'                      and optionally a \code{description} slot.
#'
#' @inherit reactome_links
setMethod("reactome_links", c("x" = "ReactomeAnalysisResult"), function(x, print_result, return_result) {
  if (missing(print_result)) print_result <- TRUE
  if (missing(return_result)) return_result <- FALSE

  if (length(x@reactome_links) == 0) {
    message("No Reactome links available\n")

    if (return_result) {
      return(list())
    }
  } else {
    # display the visualizations
    if (print_result) {
      for (i in seq(to = length(x@reactome_links))) {
        cat(i, ": ", x@reactome_links[[i]]["name"], "\n", sep = "")
        cat("   ", x@reactome_links[[i]]["description"], "\n", sep = "")
        cat("   URL = ", x@reactome_links[[i]]["url"], "\n", sep = "")
      }
    }

    # return if set
    if (return_result) {
      return(x@reactome_links)
    }
  }
})

#' open_reactome
#'
#' Opens the specified Reactome visualization in the system's default browser.
#'
#' @param x ReactomeAnalysisResult.
#' @param ... Additional parameters passed to downstream functions.
#' @return The opened link
#'
#' @export
#' @family ReactomeAnalysisResult functions
#'
#' @examples
#' # Note: This function only works with a newly created result
#' # since the visualization links only stay active for 7 days
#'
#' # load an example result
#' library(ReactomeGSA.data)
#' data(griss_melanoma_result)
#'
#' # get the reactome link - this does only work
#' # with new results
#' # open_reactome(griss_melanoma_result)
setGeneric("open_reactome", function(x, ...) standardGeneric("open_reactome"))

#' open_reactome - ReactomeAnalysisResult
#'
#' @param n_visualization numeric The index of the visualization to display (default \code{1}).
#'                        Use \code{\link{reactome_links}}
#'                        to retrieve all available visualizations and their index. By default,
#'                        the first visualization is opened.
#'
#' @inherit open_reactome
setMethod("open_reactome", c("x" = "ReactomeAnalysisResult"), function(x, n_visualization, ...) {
  if (missing(n_visualization)) n_visualization <- 1

  if (length(x@reactome_links) < 1) {
    stop("Result does not contain any visualizations.")
  }

  if (length(x@reactome_links) < n_visualization) {
    stop("Error: Visualization ", n_visualization, " not available")
  }

  url = x@reactome_links[[n_visualization]]["url"]

  message("Opening ", url, " in the default browser...")

  utils::browseURL(url = url)

  return(url)
})
