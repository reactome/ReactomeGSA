
# analyse_sc_clusters methods ----


#' analyse_sc_clusters
#' 
#' Analyses cell clusters of a single-cell RNA-sequencing
#' experiment to get pathway-level expressions for every
#' cluster of cells.
#' 
#' There are currently two specific implementations of
#' this function, one to support \code{Seurat} objects
#' and one to support Bioconductor's \code{SingleCellExperiment}
#' class.
#'
#' @param object The object containing the single-cell RNA-sequencing data.
#' @param use_interactors If set (default), protein-protein interactors from IntAct are used to
#'                            extend Reactome pathways.
#' @param include_disease_pathways If set, disease pathways are included as well. Disease pathways in
#'                                 Reactome follow a different annotation approach and can therefore
#'                                 lead to inaccurate results.
#' @param create_reactome_visualization If set, the interactive visualization in Reactome's PathwayBrowser
#'                                      is created.
#' @param create_reports If set, PDF and Microsoft Excel reports are created. Links to these report files
#'                       are send to the supplied e-mail address.
#' @param report_email The e-mail address to which reports should be sent to.
#' @param verbose If set, additional status messages are printed.
#' @param ... Parameters passed to the specific implementation. Detailed documentations
#'            can be found there.
#'
#' @return A \code{\link{ReactomeAnalysisResult}} object.
#' @export
#'
#' @examples
#' # This example shows how a Seurat object can be analysed
#' # the approach is identical for SingleCellExperiment objects
#' library(ReactomeGSA.data)
#' data(jerby_b_cells)
#' 
#' # perform the GSVA analysis
#' gsva_result <- analyse_sc_clusters(jerby_b_cells, verbose = FALSE)
setGeneric("analyse_sc_clusters", function(object, use_interactors = TRUE, 
                                           include_disease_pathways = FALSE,  
                                           create_reactome_visualization = FALSE,
                                           create_reports = FALSE,
                                           report_email = NULL,
                                           verbose = FALSE, ...) standardGeneric("analyse_sc_clusters"))

#' analyse_sc_clusters - Seurat
#'
#' @inherit analyse_sc_clusters
#' 
#' @param object The \code{Seurat} object containing the single cell RNA-sequencing data.
#' @param assay By default, the "RNA" assay is used, which contains the original read counts.
#' @param slot The slot in the Seurat object to use. Default and recommended approach is to use the raw counts.
setMethod("analyse_sc_clusters", c("object" = "Seurat"), function(object, use_interactors = TRUE, 
                                                                  include_disease_pathways = FALSE,  
                                                                  create_reactome_visualization = FALSE,
                                                                  create_reports = FALSE,
                                                                  report_email = NULL,
                                                                  verbose = FALSE, 
                                                                  assay = "RNA",
                                                                  slot = "counts", ...) {
  # make sure the assay exists
  if (!assay %in% Seurat::Assays(object)) {
    stop("Error: Assay '", assay, "' does not exist in passed Seurat object.")
  }
  
  # get the data
  raw_data <- Seurat::GetAssayData(object, assay = assay, slot = slot)
  
  # get the identis
  cell_ids <- as.character( Seurat::Idents(object) )
  
  # get the average counts
  if (verbose) message("Calculating average cluster expression...")
  
  av_counts <- apply(raw_data, 1, function(row_data) {
    by(row_data, cell_ids, mean)
  })
  
  # convert to a data.frame
  av_counts <- t( data.frame(av_counts) )
  
  # create the ReactomeGSA request
  request <- ReactomeGSA::ReactomeAnalysisRequest(method = "ssGSEA")
  
  # set the default request parameters
  request <- ReactomeGSA::set_parameters(request, use_interactors = use_interactors, 
                                         include_disease_pathways = include_disease_pathways, 
                                         create_reactome_visualization = create_reactome_visualization,
                                         create_reports = create_reports)
  
  if (!is.null(report_email)) {
    request <- ReactomeGSA::set_parameters(request, email = report_email)
  }
  
  # create the request object
  cell_groups <- colnames(av_counts)
  
  request <- ReactomeGSA::add_dataset(request, expression_values = av_counts, 
                                      name = "Seurat", type = "rnaseq_counts", 
                                      comparison_factor = "Cluster", comparison_group_1 = unique(cell_groups)[1], comparison_group_2 = unique(cell_groups)[2], 
                                      sample_data = data.frame(row.names = cell_groups, Cluster = cell_groups))
  
  # run the analysis
  gsa_res <- ReactomeGSA::perform_reactome_analysis(request, verbose = verbose)
  
  return(gsa_res)
})

#' analyse_sc_clusters - SingleCellExperiment
#'
#' @param cell_id The metadata field to use to group cells (for example "cluster")
#' @inherit analyse_sc_clusters
#' 
#' @param object The \code{SingleCellExperiment} object containing the single cell RNA-sequencing data.
#' @param ... Parameters passed to scater's \code{aggregateAcrossCells} function.
setMethod("analyse_sc_clusters", c("object" = "SingleCellExperiment"), function(object, use_interactors = TRUE, 
                                                                  include_disease_pathways = FALSE,  
                                                                  create_reactome_visualization = FALSE,
                                                                  create_reports = FALSE,
                                                                  report_email = NULL,
                                                                  verbose = FALSE,
                                                                  cell_id, ...) {
  # make sure scater is available
  if (!requireNamespace("scater")) {
    stop("Error: This function requires 'scater'. Please install it using BiocManager::install(\"scater\")")
  }
  
  # create the parameters for the AverageExpression call
  scater_params <- list(...)
  scater_params["x"] <- object
  scater_params["ids"] = cell_id
  
  # get the count data
  if (verbose) message("Calculating average expression per cluster...")
  agg_counts <- do.call(scater::aggregateAccrossCells, scater_params)
  counts <- SingleCellExperiment::counts(agg_counts)
  
  # create the ReactomeGSA request
  request <- ReactomeGSA::ReactomeAnalysisRequest(method = "ssGSEA")
  
  # set the default request parameters
  request <- ReactomeGSA::set_parameters(request, use_interactors = use_interactors, 
                                         include_disease_pathways = include_disease_pathways, 
                                         create_reactome_visualization = create_reactome_visualization,
                                         create_reports = create_reports)
  
  if (!is.null(report_email)) {
    request <- ReactomeGSA::set_parameters(request, email = report_email)
  }
  
  # create the request object
  df_counts <- data.frame(counts)
  cell_groups <- colnames(df_counts)
  
  request <- ReactomeGSA::add_dataset(request, expression_values = df_counts, 
                                      name = "Seurat", type = "rnaseq_counts", 
                                      comparison_factor = "Cluster", comparison_group_1 = unique(cell_groups)[1], comparison_group_2 = unique(cell_groups)[2], 
                                      sample_data = data.frame(row.names = cell_groups, Cluster = cell_groups))
  
  # run the analysis
  gsa_res <- ReactomeGSA::perform_reactome_analysis(request, verbose = verbose)
  
  return(gsa_res)
})