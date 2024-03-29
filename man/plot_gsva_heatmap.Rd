% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gsva_plotting_functions.R
\name{plot_gsva_heatmap}
\alias{plot_gsva_heatmap}
\title{plot_gsva_heatmap}
\usage{
plot_gsva_heatmap(
  object,
  pathway_ids = NULL,
  max_pathways = 20,
  truncate_names = TRUE,
  ...
)
}
\arguments{
\item{object}{The \code{\link{ReactomeAnalysisResult}} object.}

\item{pathway_ids}{A vector of pathway ids. If set, only these pathways are included in the plot.}

\item{max_pathways}{The maximum number of pathways to include. Only takes effect if \code{pathway_ids}
is not set.}

\item{truncate_names}{If set, long pathway names are truncated.}

\item{...}{Additional parameters passed to specific implementations.}
}
\value{
None
}
\description{
Plots pathway expression values / sample as a heatmap. Ranks pathways based on their
expression difference.
}
\examples{
# load the scRNA-seq example data
library(ReactomeGSA.data)
data(jerby_b_cells)

# perform the GSVA analysis
gsva_result <- analyse_sc_clusters(jerby_b_cells, verbose = FALSE)

# plot the heatmap
relevant_pathways <- c("R-HSA-983170", "R-HSA-388841", "R-HSA-2132295", "R-HSA-983705", "R-HSA-5690714")
plot_gsva_heatmap(gsva_result, 
                  pathway_ids = relevant_pathways, # limit to these pathways
                  margins = c(6,30), # adapt the figure margins in heatmap.2
                  dendrogram = "col", # only plot column dendrogram
                  scale = "row", # scale for each pathway
                  key = FALSE, # don't display the color key
                  lwid=c(0.1,4)) # remove the white space on the left
}
\seealso{
Other ReactomeAnalysisResult functions: 
\code{\link{get_result}()},
\code{\link{names,ReactomeAnalysisResult-method}},
\code{\link{open_reactome}()},
\code{\link{pathways}()},
\code{\link{plot_correlations}()},
\code{\link{plot_gsva_pathway}()},
\code{\link{plot_heatmap}()},
\code{\link{plot_volcano}()},
\code{\link{reactome_links}()},
\code{\link{result_types}()}
}
\concept{ReactomeAnalysisResult functions}
