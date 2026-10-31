#' Plot a clean UMAP
#'
#' Draws a UMAP coloured by a metadata column and adds compact axis arrows.
#'
#' @param seurat_obj A Seurat object.
#' @param group.by Metadata column used for colours.
#' @param reduction Reduction name. Default `"umap"`.
#' @param save_path Output PDF filename. Use `NULL` to skip writing.
#' @param palette Optional named or unnamed colour vector.
#' @param point_size Point size. Default `0.2`.
#'
#' @return A `ggplot` object.
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \dontrun{
#' p <- Plot_UMAP(seurat_obj, group.by = "cell_type")
#' }
Plot_UMAP <- function(
    seurat_obj,
    group.by = "cell_type",
    reduction = "umap",
    save_path = "UMAP_PBMC_Nature.pdf",
    palette = NULL,
    point_size = 0.2
) {
  .Mescnet_assert_seurat(seurat_obj)
  if (!group.by %in% colnames(seurat_obj@meta.data)) {
    stop("`group.by` was not found in the Seurat metadata.", call. = FALSE)
  }
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop("Reduction '", reduction, "' was not found in the Seurat object.",
         call. = FALSE)
  }

  umap_coords <- Seurat::Embeddings(seurat_obj, reduction)
  groups <- as.character(seurat_obj@meta.data[[group.by]])
  group_levels <- unique(groups)
  colors <- .Mescnet_complete_palette(palette, group_levels)
  df <- data.frame(
    UMAP_1 = umap_coords[, 1],
    UMAP_2 = umap_coords[, 2],
    cell_type = factor(groups, levels = group_levels),
    stringsAsFactors = FALSE
  )

  x_range <- range(df$UMAP_1, na.rm = TRUE)
  y_range <- range(df$UMAP_2, na.rm = TRUE)
  dx <- diff(x_range)
  dy <- diff(y_range)
  x0 <- x_range[1] + 0.07 * dx
  y0 <- y_range[1] + 0.07 * dy
  x1 <- x0 + 0.13 * dx
  y1 <- y0 + 0.13 * dy
  arrow <- grid::arrow(length = grid::unit(0.12, "cm"), type = "closed")

  p_clean <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$UMAP_1, y = .data$UMAP_2)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(color = .data$cell_type),
      size = point_size,
      alpha = 0.8
    ) +
    ggplot2::scale_color_manual(values = colors, drop = FALSE) +
    ggplot2::annotate(
      "segment", x = x0, y = y0, xend = x1, yend = y0,
      linewidth = 0.5, color = "black", arrow = arrow
    ) +
    ggplot2::annotate(
      "segment", x = x0, y = y0, xend = x0, yend = y1,
      linewidth = 0.5, color = "black", arrow = arrow
    ) +
    ggplot2::annotate(
      "text", x = x0 + 0.065 * dx, y = y0 - 0.025 * dy,
      label = "UMAP_1", size = 3.2, color = "black"
    ) +
    ggplot2::annotate(
      "text", x = x0 - 0.025 * dx, y = y0 + 0.065 * dy,
      label = "UMAP_2", angle = 90, size = 3.2, color = "black"
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 14, color = "black", face = "plain"),
      legend.position = "right",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 14, color = "black"),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(size = 4))
    )

  if (!is.null(save_path)) {
    .Mescnet_output_dir(save_path)
    ggplot2::ggsave(save_path, plot = p_clean, width = 6, height = 5)
  }
  p_clean
}
