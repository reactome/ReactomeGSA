#' Retrieves the status of the submitted analysis using \code{\link{start_reactome_analysis}}
#'
#' @param analysis_id The running analysis' id
#' @param reactome_url URL of the Reactome API Server. Overwrites the URL set in the 'reactome_gsa.url' option.
#'                     Specific ports can be set using the standard URL specification (for example http://your.service:1234)
#'
#' @return A list containing the \code{id}, \code{status} (can be "running", "complete", "failed"),
#'         \code{description}, and \code{completed} (numeric between 0 - 1)
#'
#' @export
get_reactome_analysis_status <- function(analysis_id, reactome_url) {
  if (missing(reactome_url)) reactome_url <- NULL

  reactome_url <- check_reactome_url(reactome_url)

  # get the status
  status_obj <- jsonlite::fromJSON(paste0(reactome_url, "0.1/status/", analysis_id))

  return(status_obj)
}

#' Retrieves the result of the submitted analysis using \code{\link{perform_reactome_analysis}}
#'
#' The result is only available if \code{\link{get_reactome_analysis_status}} indicates that the
#' analysis is complete.
#'
#' @param analysis_id The running analysis' id
#' @param reactome_url URL of the Reactome API Server. Overwrites the URL set in the 'reactome_gsa.url' option.
#'                     Specific ports can be set using the standard URL specification (for example http://your.service:1234)
#'
#' @return The result object
#'
#' @export
get_reactome_analysis_result <- function(analysis_id, reactome_url) {
  if (missing(reactome_url)) reactome_url <- NULL

  reactome_url <- check_reactome_url(reactome_url)

  # get the status
  result_obj <- jsonlite::fromJSON(paste0(reactome_url, "0.1/result/", analysis_id))

  # create the ReactomeAnalysisResult object
  reactome_obj <- convert_reactome_result(result_obj)

  return(reactome_obj)
}
