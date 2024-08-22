#' generate_pseudo_bulk_data 
#' @param group_by      entry in metadata table, based on these cluster annotation pseudo bulk is performed
#' @param split_by      variable -> split by a variable within the metadata; k must be a string
#'                      random -> splits based on a random number; k must be a number
#'                      Louvain, Louvain_multilevel, SLM, Leiden -> subclusters k must be a list with [resolution, cluster_1, cluster_2]
#' @param k_variable    variable dependent on the split_by
#'
#' @returns             returns pseudo bulk generated data
#' 
#' @examples
#' 
#' #using SCE object
#' library(scRNAseq)
#' SCE_OBJECT <- ZeiselBrainData()
#' # generating pseudo bulk data using the SCE object above, and clustering level level1class from the metadata
#' SCE_RESULT_RANDOM <- generate_pseudo_bulk_data(SCE_OBJECT, "level1class", "random",5)  # generate pseudo bulk data based on random subsampling
#' SCE_RESULT_VARIABLE <- generate_pseudo_bulk_data(SCE_OBJECT, "level1class","variable","tissue") # generate pseudo bulk data based on variable within the metadata
#'
#'
#'
#' @export
setGeneric("generate_pseudo_bulk_data", function(object,
                                                 group_by = NULL, 
                                                 split_by  = "random", 
                                                 k_variable = "4") standardGeneric("generate_pseudo_bulk_data"))

#' @param group_by      entry in metadata table, based on these cluster annotation pseudo bulk is performed
#' @param split_by      variable -> split by a variable within the metadata; k must be a string
#'                      random -> splits based on a random number; k must be a number
#'                      Louvain, Louvain_multilevel, SLM, Leiden -> subclusters k must be a list with [resolution, cluster_1, cluster_2]
#' @param k_variable    variable dependent on the split_by
#'
#' @returns             returns pseudo bulk generated data
#' @export
setMethod("generate_pseudo_bulk_data", c("object" = "Seurat"), function(object, 
                                                                        group_by, 
                                                                        split_by = "random", 
                                                                        k_variable){
      
  if (!is.character(split_by)) {
    stop('Error: split_by must be a string e.g "variable", "random", "Louvain",...')
  }

  if(!(group_by %in% colnames(object@meta.data))){
    stop('Error: group_by must be a column in metadata')
  }

  if(!(split_by %in% list("random","variable","Louvian","Louvain_multilevel","SLM","Leiden"))){
    stop('Error: Algorithm not found must be "random","variable","Louvian","Louvain_multilevel","SLM","Leiden"')
  }
    
  if (split_by == "variable"){
    if (!is.character(k_variable)) {
      stop('Error: k must be a string')
    }
    result <- split_variable(seurat_object, group_by, k_variable)
    return(result)
  }
  
  if (split_by == "random"){
    if (!is.numeric(k_variable)) {
      stop('Error: k must be a number')
    }
    result <- split_variable_random(seurat_object, group_by, k_variable)
    return(result)
  }
  
  if(split_by == "Louvian" || split_by== "Louvain_multilevel" || split_by == "SLM" || "Leiden"){
    if (length(k_variable) != 3) {
      stop('Error: k variables must contain [resolution, refrence cluster, comparison cluster]')
    } 

    resolution_ <- k_variable[[1]]
    subcluster_ref <- k_variable[[2]]   # subcluster variable depending on entry in metadata
    subcluster_comp <- k_variable[[3]]  # subcluster variable depending on entry in metadata


    if(!is.numeric(resolution_)){
      stop('Error: resolution must be a number')
    }

    if(split_by == "Louvain"){
      result <- split_clustering(seurat_object, group_by,resolution_, 1, subcluster_ref,subcluster_comp)
      return(result)
    }
    
    if(split_by == "Louvain_multilevel"){
      result <- split_clustering(seurat_object, group_by,resolution_, 2, subcluster_ref,subcluster_comp)
      return(result)
    }
    
    if(split_by =="SLM"){
      result <- split_clustering(seurat_object, group_by,resolution_, 3, subcluster_ref,subcluster_comp)
      return(result)
    }
    
    if(split_by =="Leiden"){
      result <- split_clustering(seurat_object, group_by,resolution_, 4, subcluster_ref,subcluster_comp)
      return(result)
    }
  }
})



#' generate_pseudo_bulk_data using SingleCellExperiment
#' @param group_by      entry in metadata table, based on these cluster annotation pseudo bulk is performed
#' @param split_by      variable -> split by a variable within the metadata; k must be a string
#'                      random -> splits based on a random number; k must be a number
#'                      subclustering [resolution, cluster_1, cluster_2]
#' @param k_variable    variable dependent on the split_by
#'
#' @returns             returns pseudo bulk generated data
#' @export
setMethod("generate_pseudo_bulk_data", c("object" = "SingleCellExperiment"), function(object, 
                                                                                      group_by, 
                                                                                      split_by, 
                                                                                      k_variable){

  if (!is.character(split_by)) {
    stop('Error: split_by must be a string e.g "variable", "random", "subclustering" ')
  }

  if(!(split_by %in% list("variable","random","subclustering"))){
    stop('Error: Algorithm not found must be "variable","random","subclustering"')
  }
    
  if (!(group_by %in% colnames(colData(object)))) {
    stop("Error: group_by must be a column in metadata")
  }

  if (split_by == "variable"){
    if (!is.character(k_variable)) {
      stop('Error: k must be a string')
    }
    result <- split_variable_sce(object, group_by, k_variable)
    return(result)
  }
  
  if (split_by == "random"){
    if (!is.numeric(k_variable)) {
      stop('Error: k must be a number')
    }
    result <- split_random_sce(object, group_by, k_variable)
    return(result)
  }
  
  if(split_by == "subclustering"){
    
    if (group_by %in% c("Louvain","Leiden","SLM","Louvain_multilevel")) {
      stop('Error: Louvain, SLM and Leiden clustering only available in Seurat')
    }
    if (length(k_variable) != 3) {
      stop('Error: k variables must contain [resolution, refrence cluster, comparison cluster]')
    } 

    resolution <- k_variable[[1]]   
    subcluster_ref <- k_variable[[2]]     # subcluster variable depending on entry in metadata
    subcluster_comp <-  k_variable[[3]]   # subcluster variable depending on entry in metadata

    if (!is.numeric(resolution)) {
      stop('Error: resolution must be a number')
    }

    result <- split_subclustering_sce(object, group_by, resolution, subcluster_ref, subcluster_comp)
    return(result)
  }
    
})




#' split SCE Object by variable
#' @param group_by      entry in metadata table, based on these cluster annotation pseudo bulk is performed
#' @param k_variable    variable for sub setting must be in the metadata
#' 
#' @returns             returns pseudo bulk generated data
#' @importFrom scuttle aggregateAcrossCells
split_variable_sce <- function(sce_object, group_by, k_variable){
  
  aggregated_object <- scuttle::aggregateAcrossCells(sce_object, ids=colData(sce_object)[,c(group_by, k_variable)])
  assay_data_aggregated <- as.data.frame(assay(aggregated_object))
  meta_data_aggregated <- colData(aggregated_object)[,c(group_by,k_variable)]
  
  clustering_level <- meta_data_aggregated[[group_by]]
  variable_pools <- meta_data_aggregated[[k_variable]]
  
  combined_list <- mapply(function(x, y) paste(x, y, sep = "_"), variable_pools, clustering_level)
  combined_list <- as.list(combined_list)
  colnames(assay_data_aggregated) <- combined_list
  
  return(assay_data_aggregated)
}



#' split SCE Object with random pooling
#' @param group_by      entry in metadata table, based on these cluster annotation pseudo bulk is performed
#' @param k_variable    number of pools that should be created
#' 
#' @returns             returns pseudo bulk generated data
#' 
#' @importFrom scuttle aggregateAcrossCells
#' 
split_random_sce <- function(sce_object, group_by, k_variable){
  metadata <- colData(sce_object)
  
  num_cells <- ncol(sce_object)
  random_data <- sample(1:k_variable, num_cells, replace = TRUE)
  colData(sce_object)$random_column <- random_data
  
  aggregated_counts <- scuttle::aggregateAcrossCells(sce_object, ids=colData(sce_object)[,c(group_by, "random_column")])
  
  meta_data_aggregated <- colData(aggregated_counts)[,c(group_by,"random_column")]
  aggregated_counts <- as.data.frame(assay(aggregated_counts))
  
  clustering_level <- meta_data_aggregated[[group_by]]
  random_pools <- meta_data_aggregated$random_column
  
  combined_list <- mapply(function(x, y) paste(x, y, sep = "_"),clustering_level, random_pools)
  combined_list <- as.list(combined_list)
  colnames(aggregated_counts) <- combined_list
  
  return(aggregated_counts)
}



#' split SCE Object with random pooling
#' @param group_by      entry in metadata table, based on these cluster annotation pseudo bulk is performed
#' 
#' @param resolution            resolution
#' @param subcluster_ref        cluster to subcluster as areference
#' @param subcluster_comp       cluster to subcluster for comparison
#' 
#' @returns             returns pseudo bulk generated data
#' 
#' @importFrom scuttle aggregateAcrossCells
#' @importFrom scran quickSubCluster
split_subclustering_sce <- function(sce_object, group_by, resolution,subcluster_ref,subcluster_comp){
  
  #check if Dim reduction and Clustering is performed 
  if(length(SingleCellExperiment::reducedDimNames(sce_object)) == 0){
    stop("No Dimensionalities for Subclustering available")
  }
  
  
  # group_by are the nn.clusters defined for subclustering -> list 
  subcluster.out <- scran::quickSubCluster(sce_object, groups=group_by,
                                    prepFUN=function(x) { # Preparing the subsetted SCE for clustering.
                                      dec <- scran::modelGeneVar(x)
                                      input <- scran::denoisePCA(x, technical=dec,
                                                          subset.row=scran::getTopHVGs(dec, prop=0.1),
                                                          BSPARAM=BiocSingular::IrlbaParam())
                                    },
                                    clusterFUN=function(x) { # Performing the subclustering in the subset.
                                      g <- scran::buildSNNGraph(x, use.dimred="PCA", k=resolution)     
                                      igraph::cluster_walktrap(g)$membership
                                    }
  )
    
  aggregated_counts_subcluster_ref <- scuttle::aggregateAcrossCells(subcluster_ref, ids=subcluster_ref$subcluster)
  aggregated_counts_subcluster_comp <- scuttle::aggregateAcrossCells(subcluster_comp, ids=subcluster_comp$subcluster)
  
  assay_data_ref <- as.data.frame(assay(aggregated_counts_subcluster_ref))
  assay_data_comp <- as.data.frame(assay(aggregated_counts_subcluster_comp))
  
  re <- cbind(assay_data_ref, assay_data_comp)
  cols <- as.list(colnames(re))  # change column names
  cols <- lapply(cols, function(x) gsub("\\.", "_", x))
  colnames(re) <- cols
  
  return (re)
}




#' split Seurat object by variable 
#' @param group_by      entry in metadata table, based on these cluster annotation pseudo bulk is performed
#' @param k_variable    variable dependent on the split_by -> meta data entry
#'
#' @returns             returns pseudo bulk generated data
#' @importFrom Seurat AggregateExpression
split_variable <- function(seurat_object, group_by, k_variable){
  seurat_object <- Seurat::AggregateExpression(seurat_object, assays = "RNA", return.seurat = T, group.by = c(group_by,k_variable))
  seurat_df <- as.data.frame(seurat_object@assays$RNA$counts)
  
    return(seurat_df) 
}


#' split Seurat object by random pooling 
#' @param group_by      entry in metadata table, based on these cluster annotation pseudo bulk is performed
#' @param k_variable    number of random pools
#'
#' @returns             returns pseudo bulk generated data
#' @importFrom Seurat AggregateExpression
split_variable_random <- function(seurat_object, group_by, k_variable){
  nrow(seurat_object@meta.data)
  seurat_object$rand_column <- sample(1:k_variable, nrow(seurat_object), replace = TRUE)
  
  seurat_object <- Seurat::AggregateExpression(seurat_object, assays = "RNA", return.seurat = T, group.by = c(group_by,"rand_column"))
  seurat_df <- as.data.frame(seurat_object@assays$RNA$counts)

  return(seurat_df)
}


#' method implementation subclustering
#' @param group_by      entry in metadata table, based on these cluster annotation pseudo bulk is performed
#' @param alg           Seurat subclustering algorithm id
#' @param cluster1      cluster to subcluster
#' @param cluster2      cluster to subcluster
#' @param k_variable    number of random pools
#'
#' @returns             returns pseudo bulk generated data
#' @importFrom Seurat FindNeighbors FindClusters AggregateExpression
split_clustering <- function(seurat_object, group_by, res, alg, cluster1, cluster2){
  cluster_ids <- list()
  cluster_ids <- append(cluster_ids, cluster1)
  cluster_ids <- append(cluster_ids, cluster2)
  
  result_dataframe <- data.frame()
  
  # setup seurat object 
  for (cluster in cluster_ids){
    seurat_object_sub = subset(seurat_object, idents = cluster, invert = FALSE)
    seurat_object_sub <- Seurat::FindNeighbors(seurat_object_sub, dim = 1:10)
    seurat_object_sub <- Seurat::FindClusters(seurat_object_sub, resolution = res, algorithm = alg, cluster.name = cluster)  
    
    # aggregate expression for cluster 
    seurat_object_sub <- Seurat::AggregateExpression(seurat_object_sub, assays = "RNA", return.seurat = T, group.by = c(group_by, cluster)) 
    
    seurat_df <- as.data.frame(seurat_object_sub@assays$RNA$counts)
    colnames(seurat_df) <- paste0(cluster, "_", 1:ncol(seurat_df)) 
    # rename columns 
    if (nrow(result_dataframe) == 0){
      result_dataframe <- data.frame(index = 1:nrow(seurat_df))
      result_dataframe <- cbind(result_dataframe,seurat_df)
    } else {
      result_dataframe <- cbind(result_dataframe,seurat_df)
    } 
  }
  
  return(result_dataframe)
}

