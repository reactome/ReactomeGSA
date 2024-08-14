#' generate_metadata 
#' @param pseudo_bulk_data  pseudobulk data generated from the generate_pseudo_bulk function 
#' @returns             meta data 
#'
#' @export
setGeneric("generate_metadata", function(pseudo_bulk_data)
  standardGeneric("generate_metadata"))
#' generate metadata
#' @param pseudo_bulk_data  pseudobulk data generated from the generate_pseudo_bulk function 
#'
#' @returns                 returns metadata table for later use
#' @export
setMethod("generate_metadata", c("pseudo_bulk_data"="data.frame"), function(pseudo_bulk_data){cols <- colnames(pseudo_bulk_data)
groups <- sapply(strsplit(cols, "_"), `[`, 1)
metadata <- data.frame(
  Group <- groups
)
colnames(metadata) <- "Group"
metadata$index_use <- colnames(pseudo_bulk_data)
rownames(metadata) <- metadata$index_use
metadata$index_use <- NULL
return(metadata)

})