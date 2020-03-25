#' Start Reactome Analysis
#'
#' Submits a \code{\link{ReactomeAnalysisRequest}} to the Reactome Analysis Service API and
#' returns the analysis id of the submitted job.
#' 
#' This function should only be used for very large requests that likely take a long time to complete.
#' By default, users should use the \code{\link{perform_reactome_analysis}} function to run an analysis.
#'
#' @param request \code{\link{ReactomeAnalysisRequest}} object to submit.
#' @param compress If set (default) the JSON request data is compressed using gzip.
#' @param reactome_url URL of the Reactome API Server. Overwrites the URL set in the 'reactome_gsa.url' option.
#'                     Specific ports can be set using the standard URL specification (for example http://your.service:1234)
#'                     
#' @export
#'                     
#' @return character The analysis job's id.
#' 
#' #' @examples
#' # create a request using Camera as an analysis
#' library(ReactomeGSA.data)
#' data(griss_melanoma_proteomics)
#'
#' my_request <- ReactomeAnalysisRequest(method = "Camera")
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
#'                           type = "proteomics_int",
#'                           comparison_factor = "condition",
#'                           comparison_group_1 = "MOCK",
#'                           comparison_group_2 = "MCM",
#'                           additional_factors = c("cell.type", "patient.id"))
#' # start the analysis      
#' analysis_id <- start_reactome_analysis(my_request)
start_reactome_analysis <- function(request, compress = TRUE, reactome_url = NULL) {

  if (!methods::is(request, "ReactomeAnalysisRequest")) {
    stop("'request' must be a ReactomeAnalysisRequest object")
  }

  # make sure the URL is OK
  reactome_url <- check_reactome_url(reactome_url = reactome_url)

  # create the json request
  json_request <- toJSON(request)
  
  # compress the JSON string
  if (compress) {
    message("Compressing request data...")
    json_request_compress <- memCompress(from = json_request, type = "gzip")
    
    # send the request
    request <- httr::POST(paste0(reactome_url, "0.1/analysis"), body = json_request_compress,
                          httr::content_type("application/gzip"), httr::add_headers("Content-Encoding" = "gzip"))
  } else {
    request <- httr::POST(paste0(reactome_url, "0.1/analysis"), body = json_request,
                          httr::content_type("application/json"))
  }

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
