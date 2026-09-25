\name{CreaterMEscnet_data}
\alias{CreaterMEscnet_data}
\title{Create a Mescnet network object from metacell counts}
\usage{
CreaterMEscnet_data(
  seurat_obj,
  nfeatures = 2000,
  group.by = "cell_type",
  reduction = "harmony",
  k = 10,
  max_shared = 5,
  target_metacells = 2000,
  group_name = "INH",
  gene_list = NULL
)
}
\arguments{
\item{seurat_obj}{A Seurat object. Must contain the dimensionality reduction
named in \code{reduction}.}

\item{nfeatures}{Number of highly variable genes to use when \code{gene_list}
is \code{NULL}. Default \code{2000}.}

\item{group.by}{Character vector of metadata columns used to define the
metacell groups. Default \code{"cell_type"}.}

\item{reduction}{Name of the dimensional reduction used to find
transcriptionally neighbouring cells. Default \code{"harmony"}.}

\item{k}{Number of nearest neighbours used to build each metacell.
Default \code{10}.}

\item{max_shared}{Maximum number of cells that two metacells may share.
Default \code{5}.}

\item{target_metacells}{Target number of metacells per group. Default
\code{2000}.}

\item{group_name}{Level(s) of \code{group.by} used to build \code{DatExpr}.
Default \code{"INH"}.}

\item{gene_list}{Optional character vector of genes to use instead of computing
highly variable genes. Default \code{NULL}.}
}
\value{
The input Seurat object with metacells, \code{DatExpr} and the network metadata
added to \code{@misc$MEscnet}.
}
\description{
\code{CreaterMEscnet_data()} mirrors \code{\link{CreaterMEscnet}} but keeps the
raw metacell count matrix instead of normalising it. Metacells are aggregated
from the \code{"data"} layer, the aggregated matrix is copied into the
\code{"data"} layer of the metacell object, and \code{DatExpr} is set from those
counts.
}
\details{
Use this variant when the downstream network should be inferred on
non-normalised metacell counts; use \code{\link{CreaterMEscnet}} for the standard
normalised workflow.
}
\seealso{
\code{\link{CreaterMEscnet}}
}
\examples{
\dontrun{
seurat_obj <- CreaterMEscnet_data(seurat_obj, group_name = "INH")
}
}
