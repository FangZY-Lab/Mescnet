#' Draw a hub-gene heat map
#'
#' Selects the top kME genes from each module and draws a Seurat heat map for
#' a subset of cell groups.
#'
#' @param seurat_obj A Seurat object.
#' @param tom_matrix Square gene-by-gene association matrix.
#' @param module_genes Named or unnamed list of gene vectors.
#' @param target_cells Metadata levels shown in the heat map.
#' @param group.by Metadata column used to select cells. Default `"cell_type"`.
#' @param top_n Number of hub genes retained per module. Default `5`.
#' @param assay Assay used for expression data. Default `"RNA"`.
#' @param palette Optional named or unnamed colour vector for cell groups.
#' @param output_file Optional output PDF. Use `NULL` to skip writing.
#'
#' @return A Seurat/ggplot heat-map object.
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \dontrun{
#' p <- Plot_Hub_Heatmap(
#'   seurat_obj,
#'   tom_matrix = tom_matrix,
#'   module_genes = module_genes,
#'   target_cells = c("INH", "EXN"),
#'   top_n = 5
#' )
#' }
Plot_Hub_Heatmap <- function(
    seurat_obj,
    tom_matrix,
    module_genes,
    target_cells,
    group.by = "cell_type",
    top_n = 5,
    assay = "RNA",
    palette = NULL,
    output_file = "Figure_E_HubGenes_Heatmap.pdf"
) {
  .Mescnet_assert_seurat(seurat_obj)
  if (!group.by %in% colnames(seurat_obj@meta.data)) {
    stop("`group.by` was not found in the Seurat metadata.", call. = FALSE)
  }
  top_n <- .Mescnet_check_top_n(top_n)
  kme_result <- .Mescnet_module_kme(
    seurat_obj, tom_matrix, module_genes, assay
  )
  module_genes <- kme_result$module_genes
  hub_genes <- .Mescnet_hub_genes(
    kme_result$kme, module_genes, kme_result$genes, top_n
  )
  if (!length(hub_genes)) {
    stop("No hub genes could be selected.", call. = FALSE)
  }

  groups <- as.character(seurat_obj@meta.data[[group.by]])
  target_cells <- as.character(target_cells)
  present_cells <- target_cells[target_cells %in% unique(groups)]
  if (!length(present_cells)) {
    stop("No cells matched `target_cells`.", call. = FALSE)
  }
  if (length(present_cells) < length(target_cells)) {
    missing_cells <- setdiff(target_cells, present_cells)
    warning(
      "These `target_cells` were not found and were skipped: ",
      paste(missing_cells, collapse = ", "),
      call. = FALSE
    )
  }
  target_cells <- present_cells
  keep_cells <- rownames(seurat_obj@meta.data)[groups %in% target_cells]
  seurat_sub <- subset(seurat_obj, cells = keep_cells)
  seurat_sub$cell_type_ordered <- factor(
    seurat_sub@meta.data[[group.by]],
    levels = target_cells
  )

  all_levels <- unique(groups)
  colors <- .Mescnet_complete_palette(palette, all_levels)
  p_heatmap <- Seurat::DoHeatmap(
    seurat_sub,
    features = hub_genes,
    assay = assay,
    disp.min = -2,
    disp.max = 2,
    angle = 45,
    group.by = "cell_type_ordered",
    group.colors = unname(colors[target_cells]),
    label = TRUE,
    draw.lines = TRUE
  ) +
    ggplot2::scale_fill_gradientn(
      colors = c("#90D7EC", "white", "#F7ACBC"),
      name = "Expression"
    ) +
    ggplot2::theme(
      text = ggplot2::element_text(size = 14, color = "black", face = "plain"),
      axis.text.y = ggplot2::element_text(
        size = 14, color = "black", face = "plain"
      ),
      legend.title = ggplot2::element_text(
        size = 14, color = "black", face = "plain"
      ),
      legend.text = ggplot2::element_text(
        size = 14, color = "black", face = "plain"
      ),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )

  if (!is.null(output_file)) {
    .Mescnet_output_dir(output_file)
    ggplot2::ggsave(output_file, plot = p_heatmap, width = 10, height = 7)
  }
  p_heatmap
}
