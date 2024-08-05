## Test code for pseudo bulk 


## Setup and Packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("scRNAseq")
BiocManager::install("scuttle")
BiocManager::install("scran")
BiocManager::install("scater")
BiocManager::install("org.Mm.eg.db", force = TRUE)


library(org.Mm.eg.db)
library(scRNAseq)
library(scuttle)
library(scran)
library(BiocSingular)
library(methods)



setGeneric("generate_pseudo_bulk_data", function(object, 
                                                 group_by = NULL, 
                                                 split_by  = "random", 
                                                 k_variable = "4") standardGeneric("generate_pseudo_bulk_data"))

#' generate_pseudo_bulk_data using Seurat object
#' @inherit generate_pseudo_bul_data 
#' @param group_by      entry in metadata table, based on these cluster annotation pseudo bulk is performed
#' @param split_by      variable -> split by a variable within the metadata; k must be a string
#'                      random -> splits based on a random number; k must be a number
#'                      Louvain, Louvain_multilevel, SLM, Leiden -> subclusters k must be a list with [resolution, cluster_1, cluster_2]
#' @param k_variable    variable dependent on the split_by
#'
#' @returns             returns pseudo bulk generated data
#' @examples
#' random <- generate_pseudo_bulk_data(test_data,"seurat_clusters","random",2)
#' variable <- generate_pseudo_bulk_data(test_data,"seurat_clusters","random","samples_id")
#' louvian_clustering <- generate_pseudo_bulk_data(test_data,"seurat_clusters","Louvian",list(4,0,1)
#'
#' @export
setMethod("generate_pseudo_bulk_data", c("object" = "Seurat"), function(object, 
                                                                        group_by, 
                                                                        split_by = "random", 
                                                                        k_variable){
  
  if (split_by == "variable"){
    #result <- split_variable(seurat_object, group_by, k_variable)
    #return(result)
    return("variable")
  }
  
  if (split_by == "random"){
    #result <- split_variable_random(seurat_object, group_by, k_variable)
    #return(result)
    return("random")
    
  }
  
  if(split_by =="Louvain"){
    resolution_ <- k_variable[[1]]
    cluster1_ <- k_variable[[2]]
    cluster2_ <- k_variable[[3]]
    
    #result <- split_clustering(seurat_object, group_by,resolution_, 1, cluster1_,cluster2_)
    #return(result)
    return("Luivan")
    
  }
  
  if(split_by == "Louvain_multilevel"){
    resolution_ <- k_variable[[1]]
    cluster1_ <- k_variable[[2]]
    cluster2_ <- k_variable[[3]]
   
     #result <- split_clustering(seurat_object, group_by,resolution_, 2, cluster1_,cluster2_)
    #return(result)
    return("laivian mulit")
    
  }
  
  if(split_by =="SLM"){
    resolution_ <- k_variable[[1]]
    cluster1_ <- k_variable[[2]]
    cluster2_ <- k_variable[[3]]
    
    #result <- split_clustering(seurat_object, group_by,resolution_, 3, cluster1_,cluster2_)
    #return(result)
    return("SLM")
    
  }
  
  if(split_by =="Leiden"){
    resolution_ <- k_variable[[1]]
    cluster1_ <- k_variable[[2]]
    cluster2_ <- k_variable[[3]]
    
    #result <- split_clustering(seurat_object, group_by,resolution_, 4, cluster1_,cluster2_)
    #return(result)
    return("Leiden")
    
    
  }
})



#' generate_pseudo_bulk_data using SingleCellExperiment
#' @inherit generate_pseudo_bulk_data 
#' @param group_by      entry in metadata table, based on these cluster annotation pseudo bulk is performed
#' @param split_by      variable -> split by a variable within the metadata; k must be a string
#'                      random -> splits based on a random number; k must be a number
#'                      subclustering [resolution, cluster_1, cluster_2]
#' @param k_variable    variable dependent on the split_by
#'
#' @returns             returns pseudo bulk generated data
#' @examples
#'
#' @export
setMethod("generate_pseudo_bulk_data", c("object" = "SingleCellExperiment"), function(object, 
                                                                                      group_by, 
                                                                                      split_by, 
                                                                                      k_variable){
  if (split_by == "variable"){
    result <- split_variable_sce(object, group_by, k_variable)
    return(result)
  }
  
  if (split_by == "random"){
    result <- split_random_sce(object, group_by, k_variable)
    return(result)
  }
  
  if(split_by == "subclustering"){ 
    result <- split_subclustering_sce(object, group_by, k_variable)
    return(result)
  }
    
})




#' split SCE Object by variable
#' @param group_by      entry in metadata table, based on these cluster annotation pseudo bulk is performed
#' @param k_variable    variable for sub setting must be in the metadata
#' 
#' @returns             returns pseudo bulk generated data
#' @examples
#' # SCE_RESULT_VARIABLE <- split_variable_sce(SCE_OBJECT, GROUP_BY,K_VARIABLE)
#' 
split_variable_sce <- function(sce_object, group_by, k_variable){
  
  aggregated_object <- aggregateAcrossCells(sce_object, ids=colData(sce_object)[,c(group_by, k_variable)])
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
#' @examples
#' # SCE_RESULT_RANDOM <- split_random_sce(SCE_OBJECT, GROUP_BY, K_NUMBER)
#' 
split_random_sce <- function(sce_object, group_by, k_variable){
  metadata <- colData(SCE_OBJECT)
  
  num_cells <- ncol(SCE_OBJECT)
  random_data <- sample(1:K_NUMBER, num_cells, replace = TRUE)
  colData(SCE_OBJECT)$random_column <- random_data
  
  aggregated_counts <- aggregateAcrossCells(SCE_OBJECT, ids=colData(SCE_OBJECT)[,c(group_by, "random_column")])
  
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
#' @param k_variable    subcluster paramters
#' 
#' @param resolution    resolution
#' @param cluster1      cluster to subcluster as areference
#' @param cluster2      cluster to subcluster for comparison
#' 
#' @returns             returns pseudo bulk generated data
#' @examples
#' K_CLUSTER <- list(3, 1,4)
#' # SCE_RESULT_CLUSTERING <- split_subclustering(SCE_OBJECT, 
#'                                                nn.clusters, 
#'                                                K_CLUSTER)
#' 
split_subclustering <- function(sce_object, group_by, k_variable){
  
  #check if Dim reduction and Clustering is performed 
  if(length(reducedDimNames(SCE_OBJECT)) == 0){
    print("No Dimensionaities for Subclustering available")
  }
  
  
  # group_by are the nn.clusters defined for subclustering -> list 
  subcluster.out <- quickSubCluster(sce_object, groups=group_by,
                                    prepFUN=function(x) { # Preparing the subsetted SCE for clustering.
                                      dec <- modelGeneVar(x)
                                      input <- denoisePCA(x, technical=dec,
                                                          subset.row=getTopHVGs(dec, prop=0.1),
                                                          BSPARAM=BiocSingular::IrlbaParam())
                                    },
                                    clusterFUN=function(x) { # Performing the subclustering in the subset.
                                      g <- buildSNNGraph(x, use.dimred="PCA", k=1)     
                                      igraph::cluster_walktrap(g)$membership  ## TODO !!!!
                                    }
  )
  
  
  subcluster_ref <- subcluster.out[[k_variable[[2]]]]
  subcluster_comp <-  subcluster.out[[k_variable[[3]]]]
  
  aggregated_counts_subcluster_ref <- aggregateAcrossCells(subcluster_ref, ids=subcluster_ref$subcluster)
  aggregated_counts_subcluster_comp <- aggregateAcrossCells(subcluster_comp, ids=subcluster_comp$subcluster)
  
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
split_variable <- function(seurat_object, group_by, k_variable){
  seurat_object <- AggregateExpression(seurat_object, assays = "RNA", return.seurat = T, group.by = c(group_by,k_variable))
  tail(Cells(seurat_object))
  seurat_df <- as.data.frame(seurat_object@assays$RNA$counts)
  
    return(seurat_df) 
}


#' split Seurat object by random pooling 
#' @param group_by      entry in metadata table, based on these cluster annotation pseudo bulk is performed
#' @param k_variable    number of random pools
#'
#' @returns             returns pseudo bulk generated data
split_variable_random <- function(seurat_object, group_by, k_variable){
  nrow(seurat_object@meta.data)
  seurat_object$rand_column <- sample(1:k_variable, nrow(seurat_object), replace = TRUE)
  
  seurat_object <- AggregateExpression(seurat_object, assays = "RNA", return.seurat = T, group.by = c(group_by,"rand_column"))
  tail(Cells(seurat_object))    
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
split_clustering <- function(seurat_object, group_by, res, alg, cluster1, cluster2){
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
  
  return(result_dataframe)
}



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

