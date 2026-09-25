\name{plot_metanet_auc}
\alias{plot_metanet_auc}
\title{Benchmark gene networks by EGAD neighbour-voting AUROC}
\usage{
plot_metanet_auc(
  result,
  cell_type = "test",
  go = .Mescnet_default_go(),
  out_dir = ".",
  width = 8,
  height = 6
)
}
\arguments{
\item{result}{A named list of square gene-gene association (or topological
overlap) matrices, with gene symbols on the row and column names. Names are used
as the network labels on the x axis, for example
\code{list(MEscnet = ..., hdwgcna = ..., random = ...)}.}

\item{cell_type}{Label used for the title of the saved figure and for the output
file name. Default \code{"test"}.}

\item{go}{Data frame of GO annotations with the GO term in column 2 and the gene
symbol in column 3. Defaults to a \code{go} object found in the global
environment (for backwards compatibility with the original scripts).}

\item{out_dir}{Directory the figure is written to. Pass \code{NULL} to skip
writing a file. Default \code{"."}.}

\item{width, height}{Figure size in inches passed to
\code{\link[ggplot2:ggsave]{ggplot2::ggsave}}. Defaults \code{8} and \code{6}.}
}
\value{
A list with

\itemize{
  \item \code{plot} - the \code{ggplot} object.
  \item \code{summary_stats} - per-network number of GO terms, minimum, maximum
    and standard deviation of the AUROC values.
  \item \code{plot_df} - long data frame of every AUROC value with its network
    label.
}

As a side effect the per-network summary statistics are printed.
}
\description{
Scores each network in \code{result} by running Gene Ontology guided gene set
analysis (GBA) with \code{\link[EGAD:run_GBA]{EGAD::run_GBA}} and plots the
resulting neighbour voting AUROC distributions as violin + notched boxplots.
}
\details{
Genes for which no Entrez ID can be mapped are dropped before the GBA run, and
the lower triangle of each matrix is restored from the upper triangle in case it
was zeroed out upstream.

Requires the suggested packages \pkg{AnnotationDbi}, \pkg{org.Hs.eg.db} and
\pkg{EGAD}.
}
\seealso{
\code{\link{run_wgcna_analysis}}
}
\examples{
\dontrun{
load("TOMs.Rdata")
res <- plot_metanet_auc(TOMs, cell_type = "INH", go = go)
res$summary_stats
}
}
