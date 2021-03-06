% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_result.R
\name{get_reactome_analysis_status}
\alias{get_reactome_analysis_status}
\title{Retrieves the status of the submitted analysis using \code{\link{start_reactome_analysis}}}
\usage{
get_reactome_analysis_status(analysis_id, reactome_url = NULL)
}
\arguments{
\item{analysis_id}{The running analysis' id}

\item{reactome_url}{URL of the Reactome API Server. Overwrites the URL set in the 'reactome_gsa.url' option.
Specific ports can be set using the standard URL specification (for example http://your.service:1234)}
}
\value{
A list containing the \code{id}, \code{status} (can be "running", "complete", "failed"),
        \code{description}, and \code{completed} (numeric between 0 - 1)
}
\description{
Retrieves the status of the submitted analysis using \code{\link{start_reactome_analysis}}
}
