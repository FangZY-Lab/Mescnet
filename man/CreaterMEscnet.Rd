\name{CreaterMEscnet}
\alias{CreaterMEscnet}
\title{Create a Mescnet network object}
\usage{
CreaterMEscnet(
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
This selects which cell population the network is inferred for. Default
\code{"INH"}.}

\item{gene_list}{Optional character vector of genes to use instead of computing
highly variable genes. Default \code{NULL}.}
}
\value{
The input Seurat object with metacells, \code{DatExpr} and the network metadata
added to \code{@misc$MEscnet}.
}
\description{
\code{CreaterMEscnet()} is the main entry point of the package. It selects the
genes used for network reconstruction, aggregates transcriptionally
neighbouring cells into metacells within each level of \code{group.by},
normalises the metacell expression matrix, and stores it as \code{DatExpr} for
the cell group of interest.
}
\details{
All results are written into the Seurat object rather than returned separately:

\itemize{
  \item \code{seurat_obj@misc$MEscnet$datExpr} - metacell expression matrix used
    for network inference.
  \item \code{seurat_obj@misc$MEscnet$MEscnet_genes} - genes selected for the
    network.
  \item \code{seurat_obj@misc$MEscnet$MEscnet_params} - parameters used for
    metacell aggregation.
  \item \code{seurat_obj@misc$MEscnet$MEscnet_metacell_obj} - the metacell Seurat
    object.
  \item \code{seurat_obj@misc$active_MEscnet} - name of the active network
    (\code{"MEscnet"}).
}

\code{CreaterMEscnet()} normalises the metacell expression matrix before
\code{DatExpr} is set. If you need the raw metacell counts instead (for example
when the metacell matrix is later used with tools that expect counts), use
\code{\link{CreaterMEscnet_data}}.

Metacell aggregation is delegated to \code{\link{MetacellsByGroups}}, the Seurat
v5 compatible replacement for \code{hdWGCNA::MetacellsByGroups()}.
}
\seealso{
\code{\link{CreaterMEscnet_data}}, \code{\link{pick_power_scale_free_power}},
\code{\link{ComputeMEscnetModules}}
}
\examples{
\dontrun{
seurat_obj <- CreaterMEscnet(
  seurat_obj,
  nfeatures = 1000,
  group.by = "All_Cells",
  reduction = "umap",
  k = 10,
  max_shared = 10,
  target_metacells = 500,
  group_name = "All_Immune",
  gene_list = common_genes
)

seurat_obj <- pick_power_scale_free_power(seurat_obj)
seurat_obj <- ComputeMEscnetModules(seurat_obj, resolution = 1, number = 1)
}
}
