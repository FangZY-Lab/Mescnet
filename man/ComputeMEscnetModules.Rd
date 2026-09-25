\name{ComputeMEscnetModules}
\alias{ComputeMEscnetModules}
\title{Detect Mescnet gene modules with Leiden community detection}
\usage{
ComputeMEscnetModules(
  seurat_obj,
  min_genes_per_module = 20,
  resolution = 1,
  number = 1,
  seed = 1234,
  objective_function = "modularity",
  n_iterations = 10
)
}
\arguments{
\item{seurat_obj}{A Seurat object containing
\code{seurat_obj@misc$MEscnet$adj_matrices}, as produced by
\code{\link{pick_power_scale_free_power}}.}

\item{min_genes_per_module}{Minimum number of genes for a community to be kept
as a module. Default \code{20}. Pass \code{0} to keep every community.}

\item{resolution}{Leiden resolution parameter; larger values yield more,
smaller modules. Default \code{1}.}

\item{number}{Index of the association matrix in
\code{seurat_obj@misc$MEscnet$adj_matrices} to use. \code{1} is the partial
correlation (\code{pcor}) network, \code{2} to \code{7} are \code{bayes_cor1} to
\code{bayes_cor6}. Default \code{1}.}

\item{seed}{Random seed used for the Leiden algorithm. Default \code{1234}.}

\item{objective_function}{Objective function passed to
\code{\link[igraph:cluster_leiden]{igraph::cluster_leiden}}. One of
\code{"modularity"} (default) or \code{"CPM"}.}

\item{n_iterations}{Number of Leiden iterations. Default \code{10}.}
}
\value{
The input Seurat object with the following added to \code{@misc$MEscnet}:

\itemize{
  \item \code{MEscnet_modules} - list with \code{module_ids},
    \code{module_sizes} and \code{gene_lists}.
  \item \code{moduleColors} - named vector over all \code{DatExpr} genes giving
    the module each gene belongs to (\code{"Module_<id>"}) or \code{"grey"}.
  \item \code{MEscnet_network} - the \code{igraph} graph that was clustered.
}
}
\description{
Turns one of the soft-thresholded adjacency matrices produced by
\code{\link{pick_power_scale_free_power}} into a weighted graph and detects gene
communities with the Leiden algorithm. Communities with at least
\code{min_genes_per_module} genes are retained as modules; the remaining genes
are labelled \code{"grey"}.
}
\details{
Leiden is used instead of Louvain-style modularity optimisation because it
guarantees well-connected communities, which avoids the poorly connected modules
that can be produced by pure modularity maximisation.
}
\seealso{
\code{\link{pick_power_scale_free_power}}, \code{\link{CreaterMEscnet}}
}
\examples{
\dontrun{
seurat_obj <- ComputeMEscnetModules(
  seurat_obj,
  min_genes_per_module = 20,
  resolution = 2,
  number = 2
)
seurat_obj@misc$MEscnet$MEscnet_modules$module_ids
}
}
