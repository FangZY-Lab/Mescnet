#' Plot the Mescnet co-expression network
#'
#' Highlights module membership and maps kME to node transparency. The top
#' hub genes in each module are labelled.
#'
#' @param seurat_obj A Seurat object.
#' @param tom_matrix Square gene-by-gene association matrix.
#' @param module_genes Named or unnamed list of gene vectors.
#' @param top_n Number of hub genes labelled per module. Default `5`.
#' @param assay Assay used for expression data. Default `"RNA"`.
#' @param edge_threshold Minimum matrix value retained as an edge.
#' @param output_file Optional output PDF. Use `NULL` to skip writing.
#'
#' @return A `ggplot` object produced with ggNetView.
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \dontrun{
#' p <- Plot_MEscnet_Network(
#'   seurat_obj,
#'   tom_matrix = tom_matrix,
#'   module_genes = module_genes,
#'   top_n = 5
#' )
#' }
Plot_MEscnet_Network <- function(
    seurat_obj,
    tom_matrix,
    module_genes,
    top_n = 5,
    assay = "RNA",
    edge_threshold = 0.02,
    output_file = "MEscnet_ggNetView_Nature_Hubs.pdf"
) {
  .Mescnet_assert_seurat(seurat_obj)
  .Mescnet_require("ggNetView")
  top_n <- .Mescnet_check_top_n(top_n)
  kme_result <- .Mescnet_module_kme(
    seurat_obj, tom_matrix, module_genes, assay
  )
  module_genes <- kme_result$module_genes
  hub_genes <- .Mescnet_hub_genes(
    kme_result$kme, module_genes, kme_result$genes, top_n
  )

  node_annotation <- data.frame(
    name = kme_result$genes,
    kME = 0,
    stringsAsFactors = FALSE
  )
  rownames(node_annotation) <- node_annotation$name
  for (i in seq_along(module_genes)) {
    genes <- intersect(module_genes[[i]], kme_result$genes)
    if (!length(genes)) {
      next
    }
    node_annotation[genes, "kME"] <- abs(
      kme_result$kme[genes, names(module_genes)[i]]
    )
  }
  node_annotation$kME[!is.finite(node_annotation$kME)] <- 0

  tom <- as.matrix(tom_matrix)
  tom[!is.finite(tom)] <- 0
  tom[tom < edge_threshold] <- 0
  diag(tom) <- 0
  tom <- (tom + t(tom)) / 2

  graph <- igraph::graph_from_adjacency_matrix(
    tom,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
  graph <- igraph::delete_vertices(
    graph,
    which(igraph::degree(graph) == 0)
  )
  if (!igraph::vcount(graph)) {
    stop("No nodes remain after applying `edge_threshold`.", call. = FALSE)
  }
  graph <- igraph::simplify(graph)
  igraph::E(graph)$weight <- abs(igraph::E(graph)$weight)
  if (is.null(igraph::V(graph)$name)) {
    igraph::V(graph)$name <- rownames(tom)
  }

  module_levels <- c(paste0("Module", seq_along(module_genes)), "Module0")
  gene_to_module <- stats::setNames(
    rep("Module0", igraph::vcount(graph)),
    igraph::V(graph)$name
  )
  for (i in seq_along(module_genes)) {
    genes <- intersect(module_genes[[i]], igraph::V(graph)$name)
    gene_to_module[genes] <- paste0("Module", i)
  }
  igraph::V(graph)$modularity2 <- unname(
    gene_to_module[igraph::V(graph)$name]
  )

  graph_table <- tidygraph::as_tbl_graph(graph)
  graph_table <- dplyr::mutate(
    graph_table,
    modularity = factor(.data$modularity2, levels = module_levels),
    modularity2 = factor(
      .data$modularity2, levels = module_levels, ordered = TRUE
    ),
    modularity3 = as.character(.data$modularity2),
    Modularity = factor(.data$modularity2, levels = module_levels),
    Degree = tidygraph::centrality_degree(mode = "out")
  )
  graph_table <- dplyr::left_join(
    graph_table, node_annotation, by = "name"
  )

  p1 <- ggNetView::ggNetView(
    graph_obj = graph_table,
    layout = "gephi",
    layout.module = "adjacent",
    group.by = "Modularity",
    fill.by = "Modularity",
    pointsize = c(2, 9),
    center = FALSE,
    jitter = FALSE,
    mapping_line = FALSE,
    shrink = 0.9,
    linealpha = 0.4,
    linecolor = "#E0E0E0"
  )

  for (i in seq_along(p1$layers)) {
    if (inherits(p1$layers[[i]]$geom, "GeomPoint")) {
      p1$layers[[i]]$aes_params$alpha <- NULL
      p1$layers[[i]]$mapping <- utils::modifyList(
        p1$layers[[i]]$mapping,
        ggplot2::aes(alpha = .data$kME)
      )
    }
  }

  module_colors <- stats::setNames(
    c(Mescnet_colors(length(module_genes)), "#CCCCCC"),
    module_levels
  )
  p_styled <- p1 +
    ggplot2::scale_fill_manual(values = module_colors, drop = FALSE) +
    ggplot2::scale_alpha_continuous(
      range = c(0.15, 1), limits = c(0, 1), name = "kME (Hubness)"
    ) +
    ggplot2::theme(
      text = ggplot2::element_text(size = 14, color = "black", face = "plain"),
      plot.title = ggplot2::element_text(
        size = 14, color = "black", face = "plain", hjust = 0.5
      ),
      legend.text = ggplot2::element_text(
        size = 14, color = "black", face = "plain"
      ),
      legend.title = ggplot2::element_text(
        size = 14, color = "black", face = "plain"
      )
    ) +
    ggplot2::labs(title = "MEscnet Co-expression Network")

  node_data <- NULL
  for (layer in p1$layers) {
    if (inherits(layer$geom, "GeomPoint")) {
      node_data <- layer$data
      break
    }
  }
  if (!is.null(node_data) && length(hub_genes)) {
    label_data <- node_data[node_data$name %in% hub_genes, , drop = FALSE]
    p_final <- p_styled +
      ggrepel::geom_text_repel(
        data = label_data,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$name),
        size = 5,
        color = "black",
        fontface = "plain",
        bg.color = "white",
        bg.r = 0.15,
        max.overlaps = Inf,
        box.padding = 0.8,
        point.padding = 0.3,
        segment.color = "grey50"
      )
  } else {
    p_final <- p_styled
  }

  if (!is.null(output_file)) {
    .Mescnet_output_dir(output_file)
    ggplot2::ggsave(output_file, plot = p_final, height = 10, width = 11)
  }
  p_final
}
