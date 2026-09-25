#' Detect Mescnet gene modules with Leiden community detection
#'
#' Turns one of the soft-thresholded adjacency matrices produced by
#' [`pick_power_scale_free_power()`] into a weighted graph and detects gene
#' communities with the Leiden algorithm. Communities with at least
#' `min_genes_per_module` genes are retained as modules; the remaining genes are
#' labelled `"grey"`.
#'
#' @param seurat_obj A Seurat object containing
#'   `seurat_obj@misc$MEscnet$adj_matrices`, as produced by
#'   [`pick_power_scale_free_power()`].
#' @param min_genes_per_module Minimum number of genes for a community to be
#'   kept as a module. Default `20`. Pass `0` to keep every community.
#' @param resolution Leiden resolution parameter; larger values yield more,
#'   smaller modules. Default `1`.
#' @param number Index of the association matrix in
#'   `seurat_obj@misc$MEscnet$adj_matrices` to use. `1` is the partial
#'   correlation (`pcor`) network, `2` to `7` are `bayes_cor1` to `bayes_cor6`.
#'   Default `1`.
#' @param seed Random seed used for the Leiden algorithm. Default `1234`.
#' @param objective_function Objective function passed to
#'   [igraph::cluster_leiden()]. One of `"modularity"` (default) or `"CPM"`.
#' @param n_iterations Number of Leiden iterations. Default `10`.
#'
#' @return The input Seurat object with the following added to
#'   `@misc$MEscnet`:
#'
#' * `MEscnet_modules` - list with `module_ids`, `module_sizes` and
#'   `gene_lists`.
#' * `moduleColors` - named vector over all `DatExpr` genes giving the module
#'   each gene belongs to (`"Module_<id>"`) or `"grey"`.
#' * `MEscnet_network` - the `igraph` graph that was clustered.
#'
#' @details
#' Leiden is used instead of Louvain-style modularity optimisation because it
#' guarantees well-connected communities, which avoids the poorly connected
#' modules that can be produced by pure modularity maximisation.
#'
#' @seealso [`pick_power_scale_free_power()`], [`CreaterMEscnet()`]
#'
#' @examples
#' \dontrun{
#' seurat_obj <- ComputeMEscnetModules(
#'   seurat_obj,
#'   min_genes_per_module = 20,
#'   resolution = 2,
#'   number = 2
#' )
#' seurat_obj@misc$MEscnet$MEscnet_modules$module_ids
#' }
#'
#' @export
ComputeMEscnetModules <- function(seurat_obj, min_genes_per_module = 20, resolution = 1,
                                  number = 1, seed = 1234,
                                  objective_function = "modularity", n_iterations = 10) {

  adj_matrices <- seurat_obj@misc$MEscnet$adj_matrices
  if (is.null(adj_matrices)) {
    stop("`seurat_obj@misc$MEscnet$adj_matrices` is NULL. Run `pick_power_scale_free_power()` first.", call. = FALSE)
  }
  if (number < 1 || number > length(adj_matrices)) {
    stop(
      "`number` must be between 1 and ", length(adj_matrices),
      " (available: ", paste(names(adj_matrices), collapse = ", "), ").",
      call. = FALSE
    )
  }
  if (!(objective_function %in% c("modularity", "CPM"))) {
    stop("`objective_function` must be either \"modularity\" or \"CPM\".", call. = FALSE)
  }

  adj <- adj_matrices[[number]]
  diag(adj) <- 0  # 去掉自环

  g <- igraph::graph_from_adjacency_matrix(adj, mode = "max", weighted = TRUE, diag = FALSE)

  igraph::V(g)$name <- rownames(adj)

  # 执行 Leiden 聚类
  set.seed(seed)
  leiden_clusters <- igraph::cluster_leiden(
    graph = g,
    objective_function = objective_function,
    resolution = resolution,
    n_iterations = n_iterations
  )

  # 获取所有节点的模块归属
  modules <- igraph::membership(leiden_clusters)

  # ---- 核心步骤：提取有效模块 ----
  # 1. 计算模块大小
  module_sizes <- table(modules)

  # 2. 筛选有效模块（>= min_genes_per_module 个基因）
  valid_module_ids <- names(module_sizes[module_sizes >= min_genes_per_module]) %>%
    as.integer()  # 注意转换类型，避免因子类型问题

  # 3. 提取模块详细信息
  valid_modules <- list(
    module_ids = valid_module_ids,
    module_sizes = module_sizes[as.character(valid_module_ids)],
    gene_lists = lapply(valid_module_ids, function(mid) {
      # 获取该模块所有基因名称（假设节点名为基因ID）
      igraph::V(g)$name[modules == mid]
    })
  )

  # 获取 datExpr 的基因名
  datExpr <- seurat_obj@misc$MEscnet$datExpr
  gene_names <- colnames(datExpr)

  # 初始化 moduleColors，长度为 datExpr 的基因数，默认 "grey"
  moduleColors <- rep("grey", length(gene_names))
  names(moduleColors) <- gene_names

  # 为每个有效模块的基因分配模块编号（Module_1, Module_2, ...）
  for (i in seq_along(valid_module_ids)) {
    module_id <- valid_module_ids[i]
    genes <- valid_modules$gene_lists[[i]]
    moduleColors[genes] <- paste0("Module_", module_id)
  }

  seurat_obj@misc$MEscnet$moduleColors <- moduleColors
  seurat_obj@misc$MEscnet$MEscnet_modules <- valid_modules
  seurat_obj@misc$MEscnet$MEscnet_network <- g

  seurat_obj
}
