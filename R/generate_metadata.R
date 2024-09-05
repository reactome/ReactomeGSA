#' generate_metadata
#' 
#' The pseudobulk data is generated using the 
#'  \code{\link{generate_pseudo_bulk()}} function.
#' 
#' @param pseudo_bulk_data Pseudobulk data generated from the `generate_pseudo_bulk` function
#' @returns Metadata table for later use
#' 
#' @examples
#' # Example pseudobulk data
#' pseudo_bulk_data <- data.frame(
#'   sample1_groupA = c(10, 20, 30),
#'   sample2_groupA = c(15, 25, 35),
#'   sample3_groupB = c(5, 10, 15)
#' )
#' 
#' # Generate metadata from pseudobulk data
#' metadata <- generate_metadata(pseudo_bulk_data) # metadata based on the group
#' 
#' @seealso \code{\link{generate_pseudo_bulk}} for generating pseudobulk data.
#' @export
setGeneric("generate_metadata", function(pseudo_bulk_data) {
  standardGeneric("generate_metadata")
})

#' Generate metadata
#' 
#' @param pseudo_bulk_data Pseudobulk data generated from the generate_pseudo_bulk function
#' 
#' @returns Returns metadata table for later use
#' 
#' @export
setMethod("generate_metadata", c("pseudo_bulk_data" = "data.frame"), function(pseudo_bulk_data) {
  cols <- colnames(pseudo_bulk_data)
  groups <- sapply(strsplit(cols, "_"), `[`, 1)
  
  metadata <- data.frame(
    Group = groups
  )
  
  colnames(metadata) <- "Group"
  metadata$index_use <- colnames(pseudo_bulk_data)
  rownames(metadata) <- metadata$index_use
  metadata$index_use <- NULL
  
  return(metadata)
})
