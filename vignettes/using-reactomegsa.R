## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
# install devtools if needed
if (!require(devtools)) {
  install.packages("devtools")
}

# install the ReactomeGSA package
if (!require(ReactomeGSA)) {
  install_github("reactome/ReactomeGSA")
}

# install the ReactomeGSA.data package
if (!require(ReactomeGSA.data)) {
  install_github("reactome/ReactomeGSA.data")
}

## ----show_methods--------------------------------------------------------
library(ReactomeGSA)

get_reactome_methods(print_methods = TRUE, return_result = FALSE)

## ----create_request------------------------------------------------------
# Create a new request object using 'Camera' for the gene set analysis
my_request <- new("ReactomeAnalysisRequest", method = "Camera")

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
                          type = "proteomics-int",
                          comparison_factor = "condition", 
                          comparison_group_1 = "MOCK", 
                          comparison_group_2 = "MCM",
                          additional_factors = c("cell.type", "patient.id"))
my_request

## ------------------------------------------------------------------------
data("griss_melanoma_rnaseq")

# only keep genes with >= 100 transcripts in total
total_transcripts <- rowSums(griss_melanoma_rnaseq$counts)
griss_melanoma_rnaseq <- griss_melanoma_rnaseq[total_transcripts >= 100, ]

# this is a edgeR DGEList object
class(griss_melanoma_rnaseq)
head(griss_melanoma_rnaseq$samples)

## ------------------------------------------------------------------------
# add the dataset
my_request <- add_dataset(request = my_request, 
                          expression_values = griss_melanoma_rnaseq, 
                          name = "RNA-seq", 
                          type = "rnaseq",
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

