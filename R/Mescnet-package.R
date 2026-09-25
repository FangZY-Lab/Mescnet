#' Mescnet: maximum-entropy single-cell gene network
#'
#' @description
#' `Mescnet` builds context-specific gene programs from single-cell and spatial
#' transcriptomics data. The workflow is:
#'
#' 1. [`CreaterMEscnet()`] aggregates transcriptionally neighbouring cells into
#'    metacells, normalises the metacell matrix and stores it as the expression
#'    matrix used for network reconstruction.
#' 2. [`pick_power_scale_free_power()`] estimates marginal and partial
#'    correlations, combines them into seven candidate `BayesCor` association
#'    matrices, and selects a soft-thresholding power for each using the
#'    scale-free topology criterion.
#' 3. [`ComputeMEscnetModules()`] turns one of those association matrices into a
#'    weighted graph and detects gene communities with the Leiden algorithm.
#' 4. [`plot_metanet_auc()`] and [`run_wgcna_analysis()`] benchmark the resulting
#'    network against reference networks using EGAD.
#'
#' The package additionally exports Seurat v5 compatible replacements for
#' `hdWGCNA::MetacellsByGroups()` and `hdWGCNA::MetaspotsByGroups()`; see
#' [`MetacellsByGroups()`] and [`MetaspotsByGroups()`].
#'
#' Results are stored in the `misc` slot of the Seurat object under
#' `seurat_obj@misc$MEscnet`, and the active network is recorded in
#' `seurat_obj@misc$active_MEscnet`.
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom magrittr %>%
NULL

# Column names used inside dplyr data-masking calls.
utils::globalVariables(c("tom", "nv_auc"))

# Abort with an actionable message when a suggested package is missing.
.Mescnet_require <- function(...) {
  pkgs <- c(...)
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "The following packages are required: ", paste(missing, collapse = ", "),
      ".\nInstall them with:\n  install.packages(",
      paste(missing, collapse = ", "), ")\n",
      "Annotation packages such as org.Hs.eg.db are installed with ",
      "BiocManager::install().",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
