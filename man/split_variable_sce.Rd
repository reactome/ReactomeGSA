% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_pseudo_bulk_data.R
\name{split_variable_sce}
\alias{split_variable_sce}
\title{split SCE Object by variable}
\usage{
split_variable_sce(sce_object, group_by, k_variable)
}
\arguments{
\item{sce_object}{The SingleCellExperiment object to analyse.}

\item{group_by}{entry in metadata table, based on these cluster annotation pseudo bulk is performed}

\item{k_variable}{variable for sub setting must be in the metadata}
}
\value{
returns pseudo bulk generated data
}
\description{
split SCE Object by variable
}
