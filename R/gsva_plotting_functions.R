#' plot_gsva_pathway
#' 
#' Plots the expression of a specific pathway from a ssGSEA result.
#'
#' @param object The \code{\link{ReactomeAnalysisResult}} object.
#' @param pathway_id The pathway's id
#' @param ... Additional parameters for specific implementations.
#'
#' @return A ggplot2 plot object
#' @export
#' @importFrom ggplot2 ggplot aes geom_point geom_hline labs
#' @family ReactomeAnalysisResult functions
#'
#' @examples
#' # load the scRNA-seq example data
#' library(ReactomeGSA.data)
#' data(jerby_b_cells)
#' 
#' # perform the GSVA analysis
#' gsva_result <- analyse_sc_clusters(jerby_b_cells, verbose = FALSE)
#' 
#' # create the plot
#' plot_obj <- plot_gsva_pathway(gsva_result, "R-HSA-389542")
setGeneric("plot_gsva_pathway", function(object, pathway_id, ...) callGeneric("plot_gsva_pathway"))

#' ReactomeAnalysisResult - plot_gsva_pathway
#'
#' @inherit plot_gsva_pathway
setMethod("plot_gsva_pathway", c("object" = "ReactomeAnalysisResult"), function(object, pathway_id, ...)  {
  # make sure it is a GSVA result
  if (!is_gsva_result(object)) {
    stop("ReactomeAnalysisResult object does not represent a GSVA result (choose method 'ssGSEA' during the analysis)")
  }
  
  gsa_pathways <- pathways(object)
  
  # get the pathway
  if (!pathway_id %in% rownames(gsa_pathways)) {
    stop("Error: Unknown pathway ", pathway_id)
  }
  
  pathway_name <- gsa_pathways[pathway_id, "Name"]
  
  # remove the name
  gsa_pathways$Name <- NULL
  
  # fix the column names
  colnames(gsa_pathways) <- gsub("\\.[^.]*$", "", colnames(gsa_pathways))
  colnames(gsa_pathways) <- gsub("\\.", " ", colnames(gsa_pathways))
  
  # get the data.frame for ggplot
  plot_data <- data.frame(
    cluster_id = colnames(gsa_pathways),
    expr = as.numeric(gsa_pathways[pathway_id, ])
  )
  
  sorted_ids <- plot_data$cluster_id[order(plot_data$expr)]
  plot_data$cluster_id <- factor(plot_data$cluster_id, levels = unique(sorted_ids))
  
  plot_obj <- ggplot2::ggplot(plot_data, ggplot2::aes(x = cluster_id, y = expr)) + 
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::labs(x = "Cluster", y = "Pathway expression", title = paste0(pathway_name, "(", pathway_id, ")")) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  return(plot_obj)
})

# plot_gsva_heatmap functions ----

#' plot_gsva_heatmap
#' 
#' Plots pathway expression values / sample as a heatmap. Ranks pathways based on their
#' expression difference.
#'
#' @param object The \code{\link{ReactomeAnalysisResult}} object.
#' @param pathway_ids A vector of pathway ids. If set, only these pathways are included in the plot.
#' @param max_pathways The maximum number of pathways to include. Only takes effect if \code{pathway_ids}
#'                     is not set.
#' @param truncate_names If set, long pathway names are truncated.
#' @param ... Additional parameters passed to specific implementations.
#'
#' @return None
#' @family ReactomeAnalysisResult functions
#' @export
#'
#' @examples
#' # load the scRNA-seq example data
#' library(ReactomeGSA.data)
#' data(jerby_b_cells)
#' 
#' # perform the GSVA analysis
#' gsva_result <- analyse_sc_clusters(jerby_b_cells, verbose = FALSE)
#' 
#' # plot the heatmap
#' relevant_pathways <- c("R-HSA-983170", "R-HSA-388841", "R-HSA-2132295", "R-HSA-983705", "R-HSA-5690714")
#' plot_gsva_heatmap(gsva_result, 
#'                   pathway_ids = relevant_pathways, # limit to these pathways
#'                   margins = c(6,30), # adapt the figure margins in heatmap.2
#'                   dendrogram = "col", # only plot column dendrogram
#'                   scale = "row", # scale for each pathway
#'                   key = FALSE, # don't display the color key
#'                   lwid=c(0.1,4)) # remove the white space on the left
setGeneric("plot_gsva_heatmap", function(object, pathway_ids = NULL, max_pathways = 20, truncate_names = TRUE, ...) callGeneric("plot_gsva_heatmap"))

#' plot_gsva_heatmap - ReactomeAnalysisResult function
#' 
#' @param ... Additional parameters passed to the heatmap.2 function.
#' 
#' @importFrom gplots heatmap.2
#' @inherit plot_gsva_heatmap
setMethod("plot_gsva_heatmap", c("object" = "ReactomeAnalysisResult"), function(object, pathway_ids = NULL, max_pathways = 20, truncate_names = TRUE, ...) {
  # make sure it is a GSVA result
  if (!is_gsva_result(object)) {
    stop("ReactomeAnalysisResult object does not represent a GSVA result (choose method 'ssGSEA' during the analysis)")
  }
  
  gsa_pathways <- pathways(object)
  
  # sort the pathways based on their difference
  max_difference <- do.call(rbind, apply(gsa_pathways, 1, function(row) {
    values <- as.numeric(row[2:length(row)])
    return(data.frame(name = row[1], min = min(values), max = max(values)))
  }))
  
  max_difference$diff <- max_difference$max - max_difference$min
  
  gsa_pathways <- gsa_pathways[order(max_difference$diff, decreasing = T), ]
  
  # get the pathways to plot
  if (!is.null(pathway_ids)) {
    # make sure all pathways exist
    if (!all(pathway_ids %in% rownames(gsa_pathways))) {
      warning("Warning: No results for the following pathways: ", paste(pathway_ids[!pathway_ids %in% rownames(gsa_pathways)], collapse = ", ", sep = ", "))
      
      pathway_ids <- pathway_ids[pathway_ids %in% rownames(gsa_pathways)]
    }
  } else {
    # use the top different pathways
    if (nrow(gsa_pathways) > max_pathways) {
      pathway_ids <- rownames(gsa_pathways)[1:max_pathways]
    } else {
      pathway_ids <- rownames(gsa_pathways)
    }
  }
  
  # get the expression values as a matrix
  expression_matrix <- as.matrix(gsa_pathways[pathway_ids, 2:ncol(gsa_pathways)])
  
  # remove the redundant dataset name
  if (length(names(object)) == 1) {
    colnames(expression_matrix) <- gsub("\\.[^\\.]*$", "", colnames(expression_matrix))
  }
  
  # use pathway names
  rownames(expression_matrix) <- gsa_pathways[rownames(expression_matrix), "Name"]
  
  # shorten the names
  if (truncate_names) {
    is_long <- nchar(rownames(expression_matrix)) > 40
    rownames(expression_matrix)[is_long] <- paste0(substr(rownames(expression_matrix)[is_long], 1, 32), "...")
  }
  
  # merge default with user parameters
  heatmap_params <- list(
    x = expression_matrix, 
    margins = c(6, 15), 
    density.info = "none", 
    trace = "none",
    col = RColorBrewer::brewer.pal(9, "RdYlBu"))
  
  user_params <- list(...)
  for (param_name in names(user_params)) {
    if (nchar(param_name) > 0) {
      heatmap_params[[param_name]] <- user_params[[param_name]]
    }
  }
  
  # create the heatmap
  do.call(gplots::heatmap.2, heatmap_params)
})

# plot_gsva_pca ----

#' plot_gsva_pca
#' 
#' Runs a Principal Component analysis (using \code{prcomp}) on the samples
#' based on the pathway analysis results.
#'
#' @param object A \code{\link{ReactomeAnalysisResult}} object containing a ssGSEA result
#' @param pathway_ids A character vector of pathway ids. If set, only these pathways will be
#'                    used for the PCA analysis.
#' @param ... Additional paramters passed to specific implementations.
#'
#' @return A ggplot2 object representing the plot.
#' @export
#'
#' @examples 
#' # load the scRNA-seq example data
#' library(ReactomeGSA.data)
#' data(jerby_b_cells)
#' 
#' # perform the GSVA analysis
#' gsva_result <- analyse_sc_clusters(jerby_b_cells, verbose = FALSE)
setGeneric("plot_gsva_pca", function(object, pathway_ids = NULL, ...) callGeneric("plot_gsva_pca"))

#' plot_gsva_pca - ReactomeAnalysisResult
#' 
#' @param ... Additional parameters are passed to \code{prcomp}
#' 
#' @inherit plot_gsva_pca
setMethod("plot_gsva_pca", c("object" = "ReactomeAnalysisResult"), function(object, pathway_ids = NULL, ...) {
  if (!is_gsva_result(object)) {
    stop("ReactomeAnalysisResult object does not represent a GSVA result (choose method 'ssGSEA' during the analysis)")
  }
  
  gsa_pathways <- pathways(object)
  
  # get the pathways to plot
  if (!is.null(pathway_ids)) {
    # make sure all pathways exist
    if (!all(pathway_ids %in% rownames(gsa_pathways))) {
      warning("Warning: No results for the following pathways: ", paste(pathway_ids[!pathway_ids %in% rownames(gsa_pathways)], collapse = ", ", sep = ", "))
      
      pathway_ids <- pathway_ids[pathway_ids %in% rownames(gsa_pathways)]
    }
  } else {
    pathway_ids <- rownames(gsa_pathways)
  }
  
  # get the expression values as a matrix
  expression_matrix <- as.matrix(gsa_pathways[pathway_ids, 2:ncol(gsa_pathways)])
  
  # remove the redundant dataset name
  if (length(names(gsva_result)) == 1) {
    colnames(expression_matrix) <- gsub("\\.[^\\.]*$", "", colnames(expression_matrix))
  }
  
  # perform the PCA
  prcomp_params <- list(x = expression_matrix)
  
  user_params = list(...)
  for (param_name in names(user_params)) {
    if (nchar(param_name) > 0 ) {
      prcomp_params[[param_name]] <- user_params[[param_name]]
    }
  }
  
  pca_fit <- do.call(stats::prcomp, prcomp_params)
  pca_summary <- summary(pca_fit)
  
  # create the data for ggplot
  plot_data <- data.frame(pca_fit$rotation)
  plot_data$sample <- rownames(pca_fit$rotation)
  
  # get the proportion of variance
  prop_var_x <- round(pca_summary$importance["Proportion of Variance", "PC1"] * 100, 1)
  prop_var_y <- round(pca_summary$importance["Proportion of Variance", "PC2"] * 100, 1)
  
  plot_obj <- ggplot2::ggplot(plot_data, ggplot2::aes(x = PC1, y = PC2, color = sample)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(x = paste0("PC1 (", prop_var_x, "% var. expl.)"), y = paste0("PC2 (", prop_var_y, "% var. expl.)")) +
    ggplot2::theme_bw()
  
  return(plot_obj)
})

#' is_gsva_result
#'
#' @param object A \code{\link{ReactomeAnalysisResult}} object
#'
#' @return Boolean indicating whether the object is a GSVA result.
#'
#' @examples
is_gsva_result <- function(object) {
  # must contain the "pathways" result
  if (!"pathways" %in% ReactomeGSA::result_types(object)) {
    return(FALSE)
  }
  
  # must not contain the FDR column
  for (dataset in names(object)) {
    if ("FDR" %in% colnames(ReactomeGSA::get_result(object, "pathways", dataset))) {
      return(FALSE)
    }
  }
  
  return(TRUE)
}