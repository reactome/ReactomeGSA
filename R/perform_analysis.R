#' Perform a REACTOME Analaysis
#'
#' This function wraps all steps required to perform
#' an Analysis using the REACTOME Analysis Service. It submits
#' the passed \code{\link{ReactomeAnalysisRequest}} object to the
#' Reactome Analysis Service API, checks the submitted analysis'
#' status and returns the result once the analysis is complete.
#'
#' @param request \code{\link{ReactomeAnalysisRequest}} to submit.
#' @param verbose logical. Indicates whether progress messages should be displayed.
#' @param reactome_url URL of the REACTOME API Server. Overwrites the URL set in the 'reactome_gsa.url' option.
#'                     Specific ports can be set using the standard URL specification (for example http://your.service:1234)
#'
#' @return The analysis' result
#' @export
#'
#' @examples
#' # create a request using Camera as an analysis
#' library(ReactomeGSA.data)
#' data(griss_melanoma_proteomics)
#'
#' my_request <- new("ReactomeAnalysisRequest", method = "Camera")
#'
#' # set maximum missing values to 0.5 and do not create any reactome visualizations
#' my_request <- set_parameters(request = my_request,
#'                              max_missing_values = 0.5,
#'                              create_reactome_visualization = FALSE)
#'
#' # add the dataset
#' my_request <- add_dataset(request = my_request,
#'                           expression_values = griss_melanoma_proteomics,
#'                           name = "Proteomics",
#'                           type = "proteomics-int",
#'                           comparison_factor = "condition",
#'                           comparison_group_1 = "MOCK",
#'                           comparison_group_2 = "MCM",
#'                           additional_factors = c("cell.type", "patient.id"))
#'
#' # perform the analysis
#' my_result <- perform_reactome_analysis(request = my_request, verbose = FALSE)
perform_reactome_analysis <- function(request, verbose = TRUE, reactome_url = NULL) {
  if (!methods::is(request, "ReactomeAnalysisRequest")) {
    stop("Error: request must be a 'ReactomeAnalysisRequest' object.")
  }

  # submit the request
  if (verbose) message("Submitting request to REACTOME API...")
  analysis_id <- start_reactome_analysis(request = request, reactome_url = reactome_url)

  # get the status
  completed <- get_reactome_analysis_status(analysis_id)

  # create a progress bar
  if (verbose)  {
    pb <- progress::progress_bar$new(total = 100, format = "Running analysis [:bar:]", show_after = 0)

    pb$message(completed[["description"]])
    pb$update(as.numeric(completed[["completed"]]))

    last_message <- completed[["description"]]
  }

  # loop until the analysis is done
  while (completed[["status"]] == "running") {
    Sys.sleep(1)
    completed <- get_reactome_analysis_status(analysis_id)

    if (verbose) {
      current_message <- completed[["description"]]

      # only update the message if it's different and the process is still running
      if (current_message != last_message && completed[["status"]] == "running") {
        pb$message(current_message)
        last_message <- current_message
      }

      pb$update(as.numeric(completed[["completed"]]))
    }
  }

  # test if the analysis failed
  if (completed[["status"]] == "failed") {
    if (verbose) warning("REACTOME Analysis failed: ", completed[["description"]])
    return(NULL)
  }

  # retrieve the result
  if (verbose) message("Retrieving result...")
  result <- get_reactome_analysis_result(analysis_id = analysis_id)

  return(result)
}
