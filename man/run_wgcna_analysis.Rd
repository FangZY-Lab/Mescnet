\name{run_wgcna_analysis}
\alias{run_wgcna_analysis}
\title{Build the Mescnet and hdWGCNA networks and benchmark them}
\usage{
run_wgcna_analysis(
  seurat_obj,
  merged_data,
  cell_type = "test",
  go = .Mescnet_default_go(),
  tom_file = "TOMs.Rdata",
  save_TOMs = TRUE
)
}
\arguments{
\item{seurat_obj}{A Seurat object produced by \code{\link{CreaterMEscnet}}, i.e.
containing \code{seurat_obj@misc$MEscnet}.}

\item{merged_data}{Metacell expression matrix with cells in rows and genes in
columns. Stored as \code{seurat_obj@misc$MEscnet$datExpr} before the analysis
runs.}

\item{cell_type}{Label used for the saved benchmark figure. Default
\code{"test"}.}

\item{go}{GO annotation table passed on to \code{\link{plot_metanet_auc}}.
Defaults to the \code{go} object in the global environment.}

\item{tom_file}{File the \code{TOMs} list is saved to. Default
\code{"TOMs.Rdata"}.}

\item{save_TOMs}{Whether to save \code{TOMs} to \code{tom_file}. Default
\code{TRUE}.}
}
\value{
A list with

\itemize{
  \item \code{seurat_obj} - the Seurat object with the hdWGCNA network added.
  \item \code{TOMs} - named list of the Mescnet, hdWGCNA and random network
    matrices.
  \item \code{plot}, \code{summary_stats}, \code{plot_df} - the benchmark outputs
    returned by \code{\link{plot_metanet_auc}}.
}
}
\description{
Runs the full network construction for a metacell expression matrix: selects
soft-thresholding powers, builds the Mescnet (partial correlation) network,
builds an hdWGCNA network from the same expression matrix, generates a random
network with the same gene set, and benchmarks all three with
\code{\link{plot_metanet_auc}}.
}
\details{
The number of metacells must exceed the number of genes, otherwise the
correlation structure is not estimable; the function stops in that case.

The genes of the hdWGCNA network are re-ordered to match the Mescnet network so
that the two matrices are directly comparable in the benchmark. The random
network is an Erdos-Renyi style symmetric matrix with the same gene set,
included as a negative control.
}
\seealso{
\code{\link{pick_power_scale_free_power}}, \code{\link{plot_metanet_auc}}
}
\examples{
\dontrun{
out <- run_wgcna_analysis(seurat_obj, merged_data, cell_type = "INH", go = go)
out$summary_stats
}
}
