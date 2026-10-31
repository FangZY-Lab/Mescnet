#' Plot UCell module scores
#'
#' Computes UCell scores for each module and draws feature plots and violin
#' plots grouped by a metadata column.
#'
#' @param seurat_obj A Seurat object.
#' @param module_genes Named or unnamed list of gene vectors.
#' @param group.by Metadata column used for the violin plots.
#' @param assay Assay used to calculate UCell scores. Default `"RNA"`.
#' @param feature_file Optional output PDF for feature plots.
#' @param violin_file Optional output PDF for violin plots.
#' @param palette Optional named or unnamed colour vector for cell groups.
#'
#' @return A list with the `Feature` and `Vln` plot objects.
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \dontrun{
#' plots <- Plot_Module_Scores(
#'   seurat_obj,
#'   module_genes = module_genes,
#'   group.by = "cell_type"
#' )
#' }
Plot_Module_Scores <- function(
    seurat_obj,
    module_genes,
    group.by = "cell_type",
    assay = "RNA",
    feature_file = "FeaturePlot_Modules.pdf",
    violin_file = "VlnPlot_Modules.pdf",
    palette = NULL
) {
  .Mescnet_assert_seurat(seurat_obj)
  .Mescnet_require("UCell")
  module_genes <- .Mescnet_module_names(module_genes)
  if (!group.by %in% colnames(seurat_obj@meta.data)) {
    stop("`group.by` was not found in the Seurat metadata.", call. = FALSE)
  }

  seurat_obj <- UCell::AddModuleScore_UCell(
    seurat_obj,
    features = module_genes,
    assay = assay
  )
  features_ucell <- paste0(names(module_genes), "_UCell")
  module_colors <- Mescnet_colors(length(features_ucell))

  feature_plots <- lapply(seq_along(features_ucell), function(i) {
    Seurat::FeaturePlot(
      seurat_obj,
      features = features_ucell[i],
      assay = assay
    ) +
      ggplot2::scale_color_gradient(
        low = "lightgrey", high = module_colors[i]
      ) +
      ggplot2::theme_void() +
      ggplot2::theme(
        text = ggplot2::element_text(
          size = 14, color = "black", face = "plain"
        ),
        plot.title = ggplot2::element_text(
          size = 14, color = "black", face = "plain", hjust = 0.5
        ),
        legend.text = ggplot2::element_text(
          size = 14, color = "black", face = "plain"
        )
      )
  })
  p_feature <- patchwork::wrap_plots(
    feature_plots,
    ncol = length(features_ucell)
  )

  group_levels <- unique(as.character(seurat_obj@meta.data[[group.by]]))
  group_colors <- .Mescnet_complete_palette(palette, group_levels)
  p_vln <- Seurat::VlnPlot(
    seurat_obj,
    features = features_ucell,
    group.by = group.by,
    assay = assay,
    cols = group_colors,
    pt.size = 0,
    ncol = 1
  ) &
    ggplot2::theme_classic() &
    ggplot2::theme(
      text = ggplot2::element_text(size = 14, color = "black", face = "plain"),
      axis.text = ggplot2::element_text(
        size = 14, color = "black", face = "plain"
      ),
      axis.title = ggplot2::element_text(
        size = 14, color = "black", face = "plain"
      ),
      plot.title = ggplot2::element_text(
        size = 14, color = "black", face = "plain", hjust = 0.5
      ),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  if (!is.null(feature_file)) {
    .Mescnet_output_dir(feature_file)
    ggplot2::ggsave(
      feature_file, plot = p_feature,
      width = 4 * length(features_ucell), height = 4
    )
  }
  if (!is.null(violin_file)) {
    .Mescnet_output_dir(violin_file)
    ggplot2::ggsave(violin_file, plot = p_vln, width = 8, height = 10)
  }

  list(Feature = p_feature, Vln = p_vln)
}
