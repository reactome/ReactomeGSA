# Functions associated with loading external datasets through ReactomeGSA ------------------------------------

#' find_public_datasets
#' 
#' Search for a public dataset in the resources supported
#' by ReactomeGSA as external data sources.
#' 
#' @param search_term The search terms as a single string. Multiple words
#'                    (seperated by a space) are combined by an "AND".
#' @param species Limit the search to selected species. The complete list of available species can be
#'                retrieved through \code{\link{get_public_species}}. By default, entries as limited to
#'                human datasets.
#' @param reactome_url URL of the Reactome API Server. Overwrites the URL set in the 'reactome_gsa.url' option.
#'                     Specific ports can be set using the standard URL specification (for example http://your.service:1234)
#' @return A data.frame containing a list of datasets found through the search.
#' @importFrom utils URLencode
#' @importFrom methods is
#' @export  
#'
#' @examples 
#' # search for any public dataset relating to BRAF in melanoma
#' melanoma_datasets <- find_public_datasets("melanoma braf")
#' 
#' # it is also possible to limit this to another species than human
#' melanoma_mouse <- find_public_datasets("melanoma", species = "Mus musculus")
#' 
#' # the list of available species can be retrieved using get_public_species
#' all_species <- get_public_species()
#' 
#' # datasets can then be loaded using the load_public_dataset function
find_public_datasets <- function(search_term, species = "Homo sapiens", reactome_url = NULL) {
    reactome_url <- check_reactome_url(reactome_url)

    # create URL safe search terms
    search_term <- utils::URLencode(search_term)
    species <- utils::URLencode(species)

    # fetch the types
    datasets <- jsonlite::fromJSON(paste0(reactome_url, "0.1/data/search?keywords=", search_term, "&species=", species))

    if(!methods::is(datasets, "data.frame")) {
        return(data.frame())
    }

    # re-order the columns
    datasets <- datasets[, c("title", "id", "resource_name", "description", "resource_loading_id", "loading_parameters", "web_link")]

    return(datasets)
}

#' load_public_dataset
#' 
#' Loads a public dataset that was found through the 
#' \code{\link{find_public_datasets}} function. The dataset
#' is returned as a Biobase ExpressionSet object.
#' 
#' @param dataset_entry The entry of the respective dataset as returned by 
#'                      the \code{\link{find_public_datasets}} function.  
#' @param verbose If set to \code{TRUE}, status messages and a status bar are displayed.
#' @param reactome_url URL of the Reactome API Server. Overwrites the URL set in the 'reactome_gsa.url' option.
#'                     Specific ports can be set using the standard URL specification (for example http://your.service:1234)
#' @return The loaded data as an ExpressionSet object.
#' @importFrom methods is
#' @export 
#' 
#' @examples 
#' # As a first step, you need to find available datasets
#' available_datasets <- find_public_datasets("psoriasis tnf")
#' 
#' # have a quick look at the found datasets
#' available_datasets[, c("id", "title")]
#' 
#' # load the first one, use the whole row of the found datasets
#' # data.frame as the parameter
#' dataset_1 <- load_public_dataset(available_datasets[1,], verbose = TRUE)
#' 
load_public_dataset <- function(dataset_entry, verbose = FALSE, reactome_url = NULL) {
    # ensure that the dataset_entry is a data.frame
    if (!methods::is(dataset_entry, "data.frame")) {
        stop("Error: 'dataset_entry' must be a row of the search result returned by 'find_public_datasets'", call. = FALSE)
    }

    # make sure the required columns are present
    if (!all(c("resource_loading_id", "loading_parameters") %in% colnames(dataset_entry)))  {
        stop("Error: 'dataset_entry' must be a row of the search result returned by 'find_public_datasets'", call. = FALSE)
    }

    # ensure that only one dataset is passed
    if (nrow(dataset_entry) != 1) {
        stop("Error: Only one dataset can be loaded at a time.", call. = FALSE)
    }

    reactome_url <- check_reactome_url(reactome_url)

    # create the URL
    loading_url <- paste0(reactome_url, "0.1/data/load/", dataset_entry[, "resource_loading_id"])

    # create the json request
    json_request <- jsonlite::toJSON(dataset_entry$loading_parameters[[1]], pretty=FALSE)
    
    # submit the request
    request <- httr::POST(loading_url, body = json_request, httr::content_type("application/json"))

    # wait until the dataset is ready
    wait_for_loading_dataset(request, verbose, reactome_url)

    # get the data
    data <- fetch_public_data(dataset_entry, reactome_url)

    return(data)
}

#' fetch_public_data
#' 
#' Loads an already available public dataset from ReactomeGSA
#' and returns it as a Biobase::ExpressionSet object.
#' 
#' @param dataset_entry The entry of the respective dataset as returned by 
#'                      the \code{\link{find_public_datasets}} function.  
#' @param reactome_url URL of the Reactome API Server. Overwrites the URL set in the 'reactome_gsa.url' option.
#'                     Specific ports can be set using the standard URL specification (for example http://your.service:1234)
#' @return The loaded data as an ExpressionSet object.
#' 
#' @import httr
#' @importClassesFrom Biobase AnnotatedDataFrame ExpressionSet MIAME
#' 
fetch_public_data <- function(dataset_entry, reactome_url) {
    # create the request
    loading_url <- paste0(reactome_url, "0.1/data/download/", dataset_entry$id)

    loaded_data <- list()

    for (data_type in c("expr", "meta")) {
        request <- httr::GET(paste0(loading_url, "?format=", data_type))

        # test if the request worked
        if (httr::status_code(request) != 200) {
            return_content <- httr::content(request, "text")
            
            if (substr(return_content, 1, 1) == "{") {
                error_obj <- jsonlite::fromJSON(return_content)
                stop("Request failed (", error_obj[["status"]], " - ", error_obj[["title"]], "): ", error_obj[["detail"]])
            } else {
                stop("Failed to download data", call. = FALSE)
            }  
        }

        # Convert the string to a data.frame
        tab_string <- httr::content(request, "text", encoding = "UTF-8")
        df <- utils::read.delim(base::textConnection(tab_string), header = TRUE)

        # save the data
        loaded_data[[data_type]] <- df
    }

    # prepare the phenoData
    sample_data <- loaded_data[["meta"]]
    rownames(sample_data) <- sample_data[, 1]
    pheno_data <- methods::new("AnnotatedDataFrame", data = sample_data)

    # prepare the expression data
    expr_data <- loaded_data[["expr"]]
    rownames(expr_data) <- expr_data[, 1]
    expr_data[, 1] <- NULL
    expr_data <- as.matrix(expr_data)

    feature_data <- methods::new("AnnotatedDataFrame", data = data.frame(GeneName = rownames(expr_data), row.names = rownames(expr_data)))

    # Create an ExpressionSet
    eset <- Biobase::ExpressionSet(assayData = expr_data, phenoData = pheno_data, featureData = feature_data)

    # Add some meta data
    miame_data <- methods::new("MIAME",
                  name = dataset_entry$id,
                  title = dataset_entry$title,
                  abstract = dataset_entry$description,
                  url = dataset_entry$web_link,
                  other = list(notes = paste0("Public dataset loaded from ", dataset_entry$resource_name, " through ReactomeGSA.")))

     Biobase::experimentData(eset) <- miame_data

    return(eset)
}

#' wait_for_loading_dataset
#' 
#' This function loops until the dataset is available. If
#' verbose is set to \code{TRUE}, the progress is displayed
#' in a status bar.
#' 
#' @param request The httr request object of the dataset loading request.
#' @param verbose If set to \code{TRUE}, the progress is displayed as a status bar.
#' @param reactome_url URL of the Reactome API Server. Overwrites the URL set in the 'reactome_gsa.url' option.
#'                     Specific ports can be set using the standard URL specification (for example http://your.service:1234)
wait_for_loading_dataset <- function(request, verbose, reactome_url) {
    # test if the request worked
    if (httr::status_code(request) != 200) {
        return_content <- httr::content(request, "text")
        
        if (substr(return_content, 1, 1) == "{") {
            error_obj <- jsonlite::fromJSON(return_content)
            stop("Request failed (", error_obj[["status"]], " - ", error_obj[["title"]], "): ", error_obj[["detail"]])
        } else {
            stop("Failed to submit dataset loading request: ", return_content)
        }  
    }

    loading_id <- httr::content(request, "text")

    # get the status
    completed <- get_dataset_loading_status(loading_id, reactome_url)

    # create a progress bar
    if (verbose)  {
        pb <- progress::progress_bar$new(total = 100, format = "Loading dataset [:bar:]", show_after = 0)

        pb$message(completed[["description"]])
        
        if (is.numeric(completed[["completed"]]))
        pb$update(as.numeric(completed[["completed"]]))

        last_message <- completed[["description"]]
    }

    # loop until the analysis is done

    # this only tracks whether the progress bar reached completion as any update afterwards
    # causes an error
    is_done <- FALSE
    error_count <- 0

    while (completed[["status"]] == "running") {
        Sys.sleep(1)
        
        completed <- tryCatch({
            get_dataset_loading_status(loading_id, reactome_url)
        },
        error=function(cond) {
            # simply ignore this the first 10 times
            if (error_count < 10) {
            error_count <- error_count + 1
            return(completed)
            }
            
            # fail if the error count is too high
            stop("Error: Failed to connect to ReactomeGSA. Please contact support if this error persists at help@reactome.org", call. = FALSE)
        })

        if (verbose) {
            current_message <- completed[["description"]]

            # only update the message if it's different and the process is still running
            if (current_message != last_message && completed[["status"]] == "running" && !is_done) {
                pb$message(current_message)
                last_message <- current_message
            }

            # only update the progress if it didn't reach "1" before, otherwise this throws an error
            if (!is_done) {
                rel_completed <- as.numeric(completed[["completed"]])

                pb$update(rel_completed)

                # prevent future progress bar updates
                if (rel_completed == 1) {
                is_done <- TRUE
                }
            }
        }
    }
}

#' get_public_species
#' 
#' Return the list of found species labels in the
#' supported public data resources
#' 
#' @param reactome_url URL of the Reactome API Server. Overwrites the URL set in the 'reactome_gsa.url' option.
#'                     Specific ports can be set using the standard URL specification (for example http://your.service:1234)
#' @return A vector of species strings.
#' 
#' @export 
#' @examples
#' 
#' # get the available species
#' available_species <- get_public_species()
#' 
#' # inspect the first 1 - 3 entries
#' available_species[1:3]
get_public_species <- function(reactome_url = NULL) {
    reactome_url <- check_reactome_url(reactome_url)

    species <- jsonlite::fromJSON(paste0(reactome_url, "0.1/data/search/species"))
    
    return(species)
}

#' Retrieves the status of the submitted dataset loading request
#'
#' @param loading_id The dataset loading process' id
#' @param reactome_url URL of the Reactome API Server. Overwrites the URL set in the 'reactome_gsa.url' option.
#'                     Specific ports can be set using the standard URL specification (for example http://your.service:1234)
#'
#' @return A list containing the \code{id}, \code{status} (can be "running", "complete", "failed"),
#'         \code{description}, and \code{completed} (numeric between 0 - 1)
#'
get_dataset_loading_status <- function(loading_id, reactome_url = NULL) {
  reactome_url <- check_reactome_url(reactome_url)

  # get the status - on error return a default status
  status_obj <- tryCatch(
    jsonlite::fromJSON(paste0(reactome_url, "0.1/data/status/", loading_id)),
    error = function(e) list(completed = 0, description = "Unknown", status = "running")
  )

  return(status_obj)
}