\name{MetacellsByGroups}
\alias{MetacellsByGroups}
\title{Aggregate single cells into metacells (Seurat v5 compatible)}
\usage{
MetacellsByGroups(
  seurat_obj,
  group.by = c("seurat_clusters"),
  ident.group = "seurat_clusters",
  k = 25,
  reduction = "pca",
  dims = NULL,
  assay = NULL,
  slot = "counts",
  layer = "counts",
  mode = "average",
  cells.use = NULL,
  min_cells = 100,
  max_shared = 15,
  target_metacells = 1000,
  max_iter = 5000,
  verbose = FALSE,
  wgcna_name = NULL
)
}
\arguments{
\item{seurat_obj}{A Seurat object.}

\item{group.by}{Character vector of metadata columns defining the metacell
groups. Default \code{"seurat_clusters"}.}

\item{ident.group}{Metadata column used to set the identity of the resulting
metacell object; must be one of \code{group.by}. Default
\code{"seurat_clusters"}.}

\item{k}{Number of nearest neighbours used to construct each metacell. Default
\code{25}.}

\item{reduction}{Name of the dimensional reduction used to find neighbouring
cells. Default \code{"pca"}.}

\item{dims}{Dimensions of \code{reduction} to use. \code{NULL} uses all
dimensions stored in the reduction.}

\item{assay}{Assay to aggregate. \code{NULL} uses the default assay.}

\item{slot, layer}{Layer or slot holding the expression values. \code{layer} is
used with SeuratObject >= 5.0.0 and \code{slot} otherwise; both default to
\code{"counts"}.}

\item{mode}{Aggregation mode, either \code{"average"} (default) or
\code{"sum"}.}

\item{cells.use}{Optional vector of cells to restrict the aggregation to.}

\item{min_cells}{Minimum number of cells a group must have to be aggregated.
Default \code{100}.}

\item{max_shared}{Maximum number of cells two metacells may share. Default
\code{15}.}

\item{target_metacells}{Target number of metacells per group. Default
\code{1000}.}

\item{max_iter}{Maximum number of iterations for the metacell construction.
Default \code{5000}.}

\item{verbose}{Print progress messages. Default \code{FALSE}.}

\item{wgcna_name}{Name of the hdWGCNA experiment the metacell object is stored
under. \code{NULL} uses \code{seurat_obj@misc$active_wgcna}.}
}
\value{
The input Seurat object with the merged metacell object and the metacell
parameters stored under \code{wgcna_name}.
}
\description{
Transcriptionally neighbouring cells within each level of \code{group.by} are
aggregated into metacells, which reduces expression sparsity and stabilises
gene-gene association estimates while retaining local cellular heterogeneity.
}
\details{
This function is a drop-in replacement for
\code{\link[hdWGCNA:MetacellsByGroups]{hdWGCNA::MetacellsByGroups()}} with two
differences:

\itemize{
  \item the expression matrix is read with the \code{layer} argument when
    \code{SeuratObject >= 5.0.0} is installed and with \code{slot} otherwise, so
    it works on both Seurat v4 and v5 objects;
  \item after merging, layers are joined with
    \code{\link[SeuratObject:JoinLayers]{SeuratObject::JoinLayers()}} on Seurat
    v5 objects, which is required before the metacell matrix can be normalised
    and used for network inference.
}

Because it shares its name with the hdWGCNA function, load \pkg{Mescnet} after
\pkg{hdWGCNA} (or call \code{Mescnet::MetacellsByGroups()}) to make sure this
version is the one that gets used.

Metacell construction relies on the internal \code{hdWGCNA:::ConstructMetacells}
function.
}
\seealso{
\code{\link{CreaterMEscnet}}, \code{\link{MetaspotsByGroups}}
}
\examples{
\dontrun{
seurat_obj <- MetacellsByGroups(
  seurat_obj,
  group.by = "cell_type",
  reduction = "harmony",
  k = 10,
  max_shared = 5,
  ident.group = "cell_type",
  target_metacells = 2000
)
}
}
