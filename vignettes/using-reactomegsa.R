## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!require(ReactomeGSA))
  BiocManager::install("ReactomeGSA")

## ----show_methods--------------------------------------------------------
library(ReactomeGSA)

available_methods <- get_reactome_methods(print_methods = FALSE, return_result = TRUE)

# only show the names of the available methods
available_methods$name

## ----get_method_details--------------------------------------------------
# Use this command to print the description of the specific method to the console
# get_reactome_methods(print_methods = TRUE, print_details = TRUE, method = "PADOG", return_result = FALSE)

# show the parameter names for the method
padog_params <- available_methods$parameters[available_methods$name == "PADOG"][[1]]

paste0(padog_params$name, " (", padog_params$type, ", ", padog_params$default, ")")

## ----create_request------------------------------------------------------
# Create a new request object using 'Camera' for the gene set analysis
my_request <-ReactomeAnalysisRequest(method = "Camera")

my_request

## ----set_parameters------------------------------------------------------
# set the maximum number of allowed missing values to 50%
my_request <- set_parameters(request = my_request, max_missing_values = 0.5)

my_request

## ----add_dataset---------------------------------------------------------
library(ReactomeGSA.data)
data("griss_melanoma_proteomics")

## ------------------------------------------------------------------------
class(griss_melanoma_proteomics)
head(griss_melanoma_proteomics$samples)

## ------------------------------------------------------------------------
my_request <- add_dataset(request = my_request, 
                          expression_values = griss_melanoma_proteomics, 
                          name = "Proteomics", 
                          type = "proteomics_int",
                          comparison_factor = "condition", 
                          comparison_group_1 = "MOCK", 
                          comparison_group_2 = "MCM",
                          additional_factors = c("cell.type", "patient.id"))
my_request

## ------------------------------------------------------------------------
data("griss_melanoma_rnaseq")

# only keep genes with >= 100 reads in total
total_reads <- rowSums(griss_melanoma_rnaseq$counts)
griss_melanoma_rnaseq <- griss_melanoma_rnaseq[total_reads >= 100, ]

# this is a edgeR DGEList object
class(griss_melanoma_rnaseq)
head(griss_melanoma_rnaseq$samples)

## ------------------------------------------------------------------------
# add the dataset
my_request <- add_dataset(request = my_request, 
                          expression_values = griss_melanoma_rnaseq, 
                          name = "RNA-seq", 
                          type = "rnaseq_counts",
                          comparison_factor = "treatment", 
                          comparison_group_1 = "MOCK", 
                          comparison_group_2 = "MCM",
                          additional_factors = c("cell_type", "patient"),
                          # This adds the dataset-level parameter 'discrete_norm_function' to the request
                          discrete_norm_function = "TMM")
my_request

## ----get_data_types------------------------------------------------------
get_reactome_data_types()

## ----perform_analysis----------------------------------------------------
result <- perform_reactome_analysis(request = my_request, verbose = TRUE)

## ------------------------------------------------------------------------
names(result)

## ------------------------------------------------------------------------
result_types(result)

## ------------------------------------------------------------------------
# retrieve the fold-change data for the proteomics dataset
proteomics_fc <- get_result(result, type = "fold_changes", name = "Proteomics")
head(proteomics_fc)

## ------------------------------------------------------------------------
combined_pathways <- pathways(result)

head(combined_pathways)

## ------------------------------------------------------------------------
sessionInfo()

