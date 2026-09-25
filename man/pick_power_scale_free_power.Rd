\name{pick_power_scale_free_power}
\alias{pick_power_scale_free_power}
\title{Select soft-thresholding powers for the Mescnet association matrices}
\usage{
pick_power_scale_free_power(
  seurat_obj,
  powerVector = 1:20,
  RsquaredCut = 0.8,
  nBreaks = 10,
  plot_diagnostics = TRUE
)
}
\arguments{
\item{seurat_obj}{A Seurat object containing
\code{seurat_obj@misc$MEscnet$datExpr}, as produced by
\code{\link{CreaterMEscnet}}.}

\item{powerVector}{Integer vector of candidate soft-thresholding powers.
Default \code{1:20}.}

\item{RsquaredCut}{R-squared threshold used to pick the recommended power.
Default \code{0.80}.}

\item{nBreaks}{Number of breaks used to bin the connectivity distribution
before fitting the log-log scale-free model. Default \code{10}.}

\item{plot_diagnostics}{Draw the four-panel diagnostic plot for each candidate
matrix. Default \code{TRUE}; set to \code{FALSE} for headless or batch runs.}
}
\value{
The input Seurat object with the following added to \code{@misc$MEscnet}:

\itemize{
  \item \code{MEscnet_ppcor} - the partial correlation matrix.
  \item \code{bayes_cor1} ... \code{bayes_cor6} - the candidate BayesCor
    association matrices.
  \item \code{result_tables} - per-matrix diagnostic table (R-squared, slope,
    mean, median and maximum connectivity for every power).
  \item \code{recommended_powers} - named list of the selected power per matrix,
    \code{NA} when no power reached \code{RsquaredCut}.
  \item \code{adj_matrices} - named list of soft-thresholded adjacency matrices
    for the matrices with a selected power.
  \item \code{tom} - the adjacency matrix of the \code{pcor} network, used by
    \code{\link{ComputeMEscnetModules}} as the default network.
}

When \code{plot_diagnostics = TRUE}, a four-panel diagnostic plot (R-squared,
mean, median and maximum connectivity versus power) is drawn for each candidate
matrix as a side effect.
}
\description{
Builds seven candidate gene-gene association matrices from the metacell
expression matrix (\code{seurat_obj@misc$MEscnet$datExpr}), evaluates the
scale-free topology fit of each over \code{powerVector}, and records the first
power whose R-squared exceeds \code{RsquaredCut}.
}
\details{
The seven matrices combine the two statistical views of gene-gene relationships
that motivate Mescnet: the marginal Pearson correlation (\code{cor}) and the
partial correlation (\code{pcor}), which estimates conditional associations
given all other genes.

The six BayesCor scores are defined as:

\itemize{
  \item \code{bayes_cor1}: \code{abs(pcor) * abs(cor)}.
  \item \code{bayes_cor2}: \code{pmax(pcor, 0) * pmax(cor, 0)}.
  \item \code{bayes_cor3}: min-max normalised \code{pcor * cor}.
  \item \code{bayes_cor4}: signed-hybrid normalised \code{pcor * cor}, i.e.
    rescaling both matrices to \code{[0, 1]} with \code{(1 + x) / 2}.
  \item \code{bayes_cor5}: z-scored \code{pcor * cor}.
  \item \code{bayes_cor6}: \code{((1 + pmax(pcor, 0)) / 2) * pmax(cor, 0)}.
}

\code{NA} values produced by zero-variance genes are replaced by \code{0}.
}
\seealso{
\code{\link{CreaterMEscnet}}, \code{\link{ComputeMEscnetModules}}
}
\examples{
\dontrun{
seurat_obj <- pick_power_scale_free_power(seurat_obj)
seurat_obj@misc$MEscnet$recommended_powers
}
}
