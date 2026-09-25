\name{MetaspotsByGroups}
\alias{MetaspotsByGroups}
\title{Aggregate spatial spots into metaspots (Seurat v5 compatible)}
\usage{
MetaspotsByGroups(
  seurat_obj,
  group.by = c("seurat_clusters"),
  ident.group = "seurat_clusters",
  assay = "Spatial",
  slot = "counts",
  layer = "counts",
  mode = "sum",
  min_spots = 50,
  wgcna_name = NULL
)
}
\arguments{
\item{seurat_obj}{A Seurat object with spatial coordinates stored in
\code{seurat_obj@meta.data} as \code{row}, \code{col}, \code{imagerow} and
\code{imagecol}.}

\item{group.by}{Character vector of metadata columns defining the groups to
aggregate within. Default \code{"seurat_clusters"}.}

\item{ident.group}{Metadata column used to set the identity of the resulting
metaspot object; must be one of \code{group.by}. Default
\code{"seurat_clusters"}.}

\item{assay}{Spatial assay to aggregate. Default \code{"Spatial"}.}

\item{slot, layer}{Layer or slot holding the expression values. \code{layer} is
used with Seurat v5 and \code{slot} otherwise; both default to
\code{"counts"}.}

\item{mode}{Aggregation mode, either \code{"sum"} (default) or
\code{"average"}.}

\item{min_spots}{Minimum number of spots a group must have to be aggregated.
Default \code{50}.}

\item{wgcna_name}{Name of the hdWGCNA experiment the metaspot object is stored
under. \code{NULL} uses \code{seurat_obj@misc$active_wgcna}.}
}
\value{
The input Seurat object with the merged metaspot object stored under
\code{wgcna_name}.
}
\description{
Merges spatially adjacent spots within each level of \code{group.by} into
metaspots, the spatial analogue of metacells, to reduce sparsity in spatial
transcriptomics data.
}
\details{
This function is a drop-in replacement for
\code{\link[hdWGCNA:MetaspotsByGroups]{hdWGCNA::MetaspotsByGroups()}}. Compared
with the hdWGCNA version it additionally:

\itemize{
  \item writes the \code{ident.group} identity of each source group back into
    the metadata of the corresponding metaspot object, which is required for
    \code{Idents<-} to succeed on Seurat v5 objects;
  \item keeps \code{groupings} and the metaspot list in sync when a group fails
    to aggregate, instead of silently dropping only one of the two;
  \item joins assay layers after merging on Seurat v5 objects.
}

Metaspot construction relies on the internal
\code{hdWGCNA:::ConstructMetaspots} function.
}
\seealso{
\code{\link{MetacellsByGroups}}
}
\examples{
\dontrun{
seurat_obj <- MetaspotsByGroups(
  seurat_obj,
  group.by = "cell_type",
  ident.group = "cell_type",
  assay = "Spatial",
  min_spots = 50
)
}
}
