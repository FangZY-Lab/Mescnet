\name{Mescnet-package}
\alias{Mescnet}
\alias{Mescnet-package}
\docType{package}
\title{Mescnet: maximum-entropy single-cell gene network}
\description{
\code{Mescnet} builds context-specific gene programs from single-cell and
spatial transcriptomics data. The workflow is:

\enumerate{
  \item \code{\link{CreaterMEscnet}} aggregates transcriptionally neighbouring
    cells into metacells, normalises the metacell matrix and stores it as the
    expression matrix used for network reconstruction.
  \item \code{\link{pick_power_scale_free_power}} estimates marginal and partial
    correlations, combines them into seven candidate \code{BayesCor} association
    matrices, and selects a soft-thresholding power for each using the
    scale-free topology criterion.
  \item \code{\link{ComputeMEscnetModules}} turns one of those association
    matrices into a weighted graph and detects gene communities with the Leiden
    algorithm.
  \item \code{\link{plot_metanet_auc}} and \code{\link{run_wgcna_analysis}}
    benchmark the resulting network against reference networks using EGAD.
}

The package additionally exports Seurat v5 compatible replacements for
\code{hdWGCNA::MetacellsByGroups()} and \code{hdWGCNA::MetaspotsByGroups()};
see \code{\link{MetacellsByGroups}} and \code{\link{MetaspotsByGroups}}.

Results are stored in the \code{misc} slot of the Seurat object under
\code{seurat_obj@misc$MEscnet}, and the active network is recorded in
\code{seurat_obj@misc$active_MEscnet}.
}
\details{
The typical analysis is:

\preformatted{
seurat_obj <- CreaterMEscnet(seurat_obj, gene_list = genes, group_name = "INH")
seurat_obj <- pick_power_scale_free_power(seurat_obj)
seurat_obj <- ComputeMEscnetModules(seurat_obj, resolution = 1, number = 1)

modules <- seurat_obj@misc$MEscnet$MEscnet_modules
colors  <- seurat_obj@misc$MEscnet$moduleColors
datExpr <- seurat_obj@misc$MEscnet$datExpr
}

\code{number = 1} selects the partial correlation network; \code{2} to
\code{7} select \code{bayes_cor1} to \code{bayes_cor6}.
}
\seealso{
\code{\link{CreaterMEscnet}}, \code{\link{CreaterMEscnet_data}},
\code{\link{pick_power_scale_free_power}},
\code{\link{ComputeMEscnetModules}}, \code{\link{plot_metanet_auc}},
\code{\link{run_wgcna_analysis}}, \code{\link{MetacellsByGroups}},
\code{\link{MetaspotsByGroups}}
}
\keyword{internal}
