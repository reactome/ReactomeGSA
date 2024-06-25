#' analyse pseudo bulk
#' In the follwoing R script the methods for generating pseudo bulk data from scRNA-seq data is provided
#' It includes the supsamping method to generate pseudobulk data later used by existing methods to perform pathway analysis




#' method implementation  variable
#' @param group_by      entry in metadata table, based on these cluster annotation pseudo bulk is performed
#' @param k_variable    variable dependant on the split_by -> meta data entry
#'
#' @returns             returns pseudo bulk generated data
#' @export
split_variable <- function(seurat_object, group_by, k_variable){
    log_info("Splitting Variable")

    seurat_object <- AggregateExpression(seurat_object, assays = "RNA", return.seurat = T, group.by = c(group_by,k_variable))
    tail(Cells(seurat_object))
    seurat_df <- as.data.frame(seurat_object@assays$RNA$counts)

    log_info('Retrive Result')
    return(seurat_df) 
}


#' method implementation random
#' @param group_by      entry in metadata table, based on these cluster annotation pseudo bulk is performed
#' @param k_variable    number of random pools
#'
#' @returns             returns pseudo bulk generated data
#' @export
split_variable_random <- function(seurat_object, group_by, k_variable){
    log_info("Splitting Random")
    
    nrow(seurat_object@meta.data)
    seurat_object$rand_column <- sample(1:k_variable, nrow(seurat_object), replace = TRUE)

    seurat_object <- AggregateExpression(seurat_object, assays = "RNA", return.seurat = T, group.by = c(group_by,"rand_column"))
    tail(Cells(seurat_object))    
    seurat_df <- as.data.frame(seurat_object@assays$RNA$counts)
        log_info('Retrive Result')

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
#' @export
split_clustering <- function(seurat_object, group_by, res, alg, cluster1, cluster2){
    log_info("SubClustering")

    cluster_ids <- list()
    cluster_ids <- append(cluster_ids, cluster1)
    cluster_ids <- append(cluster_ids, cluster2)
        
    result_dataframe <- data.frame()

    # setup seurat object 
    Idents(seurat_object) <- group_by     
    for (cluster in cluster_ids){
        seurat_object_sub = subset(seurat_object, idents = cluster, invert = FALSE) 
        seurat_object_sub <- FindNeighbors(seurat_object_sub, dim = 1:10)
        seurat_object_sub <- FindClusters(seurat_object_sub, resolution = res, algorithm = alg, cluster.name = cluster)  

        # aggregate expression for cluster 
        seurat_object_sub <- AggregateExpression(seurat_object_sub, assays = "RNA", return.seurat = T, group.by = c(group_by, cluster)) 

        tail(Cells(seurat_object))    
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
    
    log_info('Retrive Result')

    return(result_dataframe)
}


#' generate_pseudo_bulk_data
#' @param group_by      entry in metadata table, based on these cluster annotation pseudo bulk is performed
#' @param split_by      variable -> split by a variable wihtin the metadata; k must be a string
#'                      random -> splits based on a random number; k must be a number
#'                      Louvian, Louvian_multilevel, SLM, Leiden -> subclusters k must be a list with [resolution, cluster_1, cluster_2]
#' @param k_variable    variable dependant on the split_by
#'
#' @returns             returns pseudo bulk generated data
#' @examples
#' random <- generate_pseudo_bulk_data(test_data,"seurat_clusters","random",2)
#' variable <- generate_pseudo_bulk_data(test_data,"seurat_clusters","random","samples_id")
#' louvian_clustering <- generate_pseudo_bulk_data(test_data,"seurat_clusters","Louvian",list(4,0,1)
#'
#' @export
generate_pseudo_bulk_data <- function(seurat_object, group_by, split_by = "random", k_variable = "3") {
    log_threshold(INFO)

    log_info('Perform Analysis:')
    log_info('Group by: {group_by}')
    log_info('Split by: {split_by}')
    log_info('k: {k_variable}')

    
    if (split_by == "variable"){
        log_info("Split by variable")
        result <- split_variable(seurat_object, group_by, k_variable)
        return(result)
    }

    if (split_by == "random"){
        log_info("Split random")
        result <- split_variable_random(seurat_object, group_by, k_variable)
        return(result)
    }

    if(split_by =="Louvian"){
        log_info("Subclustering Louvian")
        resolution_ <- k_variable[[1]]
        cluster1_ <- k_variable[[2]]
        cluster2_ <- k_variable[[3]]
        
        result <- split_clustering(seurat_object, group_by,resolution_, 1, cluster1_,cluster2_)
        return(result)
    }

    if(split_by == "Louvian_multilevel"){
        log_info("Subclustering Louvian Mulitlevel")
        resolution_ <- k_variable[[1]]
        cluster1_ <- k_variable[[2]]
        cluster2_ <- k_variable[[3]]
        result <- split_clustering(seurat_object, group_by,resolution_, 2, cluster1_,cluster2_)
        return(result)
    }

    if(split_by =="SLM"){
        log_info("Subclustering  SLM")
        resolution_ <- k_variable[[1]]
        cluster1_ <- k_variable[[2]]
        cluster2_ <- k_variable[[3]]
        
        result <- split_clustering(seurat_object, group_by,resolution_, 3, cluster1_,cluster2_)
        return(result)
    }

    if(split_by =="Leiden"){
        log_info("Subclustering  Leiden")
        resolution_ <- k_variable[[1]]
        cluster1_ <- k_variable[[2]]
        cluster2_ <- k_variable[[3]]
        
        result <- split_clustering(seurat_object, group_by,resolution_, 4, cluster1_,cluster2_)
        return(result)

    }
}

### function to automatically generate the metadata ###

#' generate metadata
#' @param pseudo_bulk_data  pseudobulk data generated from the generate_pseudo_bulk function 
#'
#' @returns                 returns metadata table for later use
#' @export
generate_metadata <- function(pseudo_bulk_data) {
    cols <- colnames(pseudo_bulk_data)
    groups <- sapply(strsplit(cols, "_"), `[`, 1)
    metadata <- data.frame(
        Group <- groups
    )
    colnames(metadata) <- "Group"
    
    metadata$index_use <- colnames(pseudo_bulk_data)
    rownames(metadata) <- metadata$index_use
    metadata$index_use <- NULL
    return(metadata)
}