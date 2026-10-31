# Mescnet 0.2.0

## New plotting API

* Added `Plot_UMAP()` for publication-style UMAPs with axis arrows and optional
  named palettes.
* Added `Plot_Module_Scores()` for UCell module scores, feature plots, and
  grouped violin plots.
* Added `Plot_MEscnet_Network()` for module-coloured networks with kME-scaled
  node transparency and per-module hub-gene labels.
* Added `Plot_Hub_Heatmap()` for heat maps of the top hub genes from each
  module across selected cell groups.
* Added `Mescnet_colors()` for a consistent high-contrast plotting palette.

## Case study and documentation

* Added a runnable PBMC case study in `case-study/`, including the master
  script, parameter guide, input-data instructions, and all result figures.
* Reworked the README so installation is the first step and dependencies are
  organized by task. hdWGCNA is documented as one internal dependency rather
  than the organizing principle of the package.
* Added an explicit mapping from `ComputeMEscnetModules(number = )` to the
  association matrix used for module detection.
* Clarified that `top_n` is the number of hub genes selected per module.

## Fixes

* Corrected the PBMC example so `module_genes1` is defined from
  `seurat_obj@misc$MEscnet$MEscnet_modules$gene_lists` before plotting.

# Mescnet 0.1.0

First release, packaging the network-construction and benchmarking functions
used in the Mescnet manuscript.

## Public API

* `CreaterMEscnet()` and `CreaterMEscnet_data()` build the metacell object and
  set `DatExpr`.
* `pick_power_scale_free_power()` builds the partial-correlation network and six
  `BayesCor` association matrices and selects soft-thresholding powers.
* `ComputeMEscnetModules()` detects modules with Leiden community detection.
* `plot_metanet_auc()` and `run_wgcna_analysis()` benchmark networks with EGAD.
* `MetacellsByGroups()` and `MetaspotsByGroups()` are Seurat v5 compatible
  replacements for the hdWGCNA functions of the same name.

## Fixes relative to the original scripts

* `run_wgcna_analysis()`: the Mescnet network is now placed in `TOMs[[1]]`
  before the hdWGCNA matrix is relabelled. Previously the code read
  `colnames(TOMs[[1]])` from an empty list, which aborted the run.
* `pick_power_scale_free_power()`: diagnostic plots can be switched off with
  `plot_diagnostics = FALSE`, and a missing `DatExpr` now raises an actionable
  error instead of a `NULL` subscript error.
* `ComputeMEscnetModules()`: `module_sizes` is now looked up by module id rather
  than by position, keeping sizes aligned with `module_ids`. The random seed,
  the Leiden objective function and the iteration count are exposed as
  arguments, `number` is checked against the available networks, and a missing
  `adj_matrices` entry raises an actionable error.
* `plot_metanet_auc()`: the GO annotation table is an explicit `go` argument
  that falls back to the global `go` object; the figure directory and size are
  configurable; `ggplot2` line width arguments are used instead of the
  deprecated `size`.
* `MetacellsByGroups()`: the `max_shared` warning now reports the user-supplied
  value instead of the value after it has been reset, groups that fail to
  aggregate are dropped, and the function stops with a clear message when every
  group fails rather than erroring later inside `rbind()`.
* `MetaspotsByGroups()`: `groupings` and the metaspot list are kept in sync when
  a group fails, layers are joined on Seurat v5, and the `ident.group` identity
  is written back into the metaspot metadata.

## Dependencies

* Core network construction needs `Seurat`, `SeuratObject`, `hdWGCNA`,
  `igraph`, `ppcor`, `ggplot2`, `dplyr`, `magrittr` and `viridis`.
* `EGAD`, `AnnotationDbi` and `org.Hs.eg.db` are only needed for
  `plot_metanet_auc()` and `run_wgcna_analysis()` and are therefore suggested
  rather than imported.
