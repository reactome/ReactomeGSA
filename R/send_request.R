#' Start Reactome Analysis
#'
#' Submits a \code{\link{ReactomeAnalysisRequest}} to the Reactome Analysis Service API and
#' returns the analysis id of the submitted job.
#'
#' @param request \code{\link{ReactomeAnalysisRequest}} object to submit.
#' @param reactome_url URL of the Reactome API Server. Overwrites the URL set in the 'reactome_gsa.url' option.
#'                     Specific ports can be set using the standard URL specification (for example http://your.service:1234)
#' @return character The analysis job's id.
start_reactome_analysis <- function(request, reactome_url = NULL) {

  if (!methods::is(request, "ReactomeAnalysisRequest")) {
    stop("'request' must be a ReactomeAnalysisRequest object")
  }

  # make sure the URL is OK
  reactome_url <- check_reactome_url(reactome_url = reactome_url)

  # create the json request
  json_request <- toJSON(request)
  
  # compress the JSON string
  json_request_compress <- memCompress(from = json_request, type = "gzip")

  # send the request
  request <- httr::POST(paste0(reactome_url, "0.1/analysis"), body = json_request_compress,
                        httr::content_type("application/gzip"), httr::add_headers("Content-Encoding" = "gzip"))

  # test if the request worked
  if (httr::status_code(request) != 200) {
    return_content <- httr::content(request, "text")
    
    if (substr(return_content, 1, 1) == "{") {
      error_obj <- jsonlite::fromJSON(return_content)
      stop("Request failed (", error_obj[["status"]], " - ", error_obj[["title"]], "): ", error_obj[["detail"]])
    } else {
      stop("Failed to submit analysis request: ", return_content)
    }  
  }

  analysis_id = httr::content(request, "text")

  message("Reactome Analysis submitted succesfully")

  # return the analysis_id
  return(analysis_id)
}
