#' Build the Mescnet and hdWGCNA networks and benchmark them
#'
#' Runs the full network construction for a metacell expression matrix: selects
#' soft-thresholding powers, builds the Mescnet (partial correlation) network,
#' builds an [hdWGCNA] network from the same expression matrix, generates a
#' random network with the same gene set, and benchmarks all three with
#' [`plot_metanet_auc()`].
#'
#' @param seurat_obj A Seurat object produced by [`CreaterMEscnet()`], i.e.
#'   containing `seurat_obj@misc$MEscnet`.
#' @param merged_data Metacell expression matrix with cells in rows and genes in
#'   columns. Stored as `seurat_obj@misc$MEscnet$datExpr` before the analysis
#'   runs.
#' @param cell_type Label used for the saved benchmark figure. Default
#'   `"test"`.
#' @param go GO annotation table passed on to [`plot_metanet_auc()`]. Defaults to
#'   the `go` object in the global environment.
#' @param tom_file File the `TOMs` list is saved to. Default `"TOMs.Rdata"`.
#' @param save_TOMs Whether to save `TOMs` to `tom_file`. Default `TRUE`.
#'
#' @return A list with
#'
#' * `seurat_obj` - the Seurat object with the hdWGCNA network added.
#' * `TOMs` - named list of the Mescnet, hdWGCNA and random network matrices.
#' * `plot`, `summary_stats`, `plot_df` - the benchmark outputs returned by
#'   [`plot_metanet_auc()`].
#'
#' @details
#' The number of metacells must exceed the number of genes, otherwise the
#' correlation structure is not estimable; the function stops in that case.
#'
#' The genes of the hdWGCNA network are re-ordered to match the Mescnet network
#' so that the two matrices are directly comparable in the benchmark. The
#' random network is an Erdos-Renyi style symmetric matrix with the same gene
#' set, included as a negative control.
#'
#' @seealso [`pick_power_scale_free_power()`], [`plot_metanet_auc()`]
#'
#' @examples
#' \dontrun{
#' out <- run_wgcna_analysis(seurat_obj, merged_data, cell_type = "INH", go = go)
#' out$summary_stats
#' }
#'
#' @export
run_wgcna_analysis <- function(seurat_obj, merged_data, cell_type = "test", go = .Mescnet_default_go(),
                               tom_file = "TOMs.Rdata", save_TOMs = TRUE) {

  # 检查维度
  dim_info <- dim(merged_data)
  if (is.null(dim_info) || length(dim_info) != 2) {
    stop("`merged_data` must be a two-dimensional matrix with cells in rows and genes in columns.", call. = FALSE)
  }
  sample_count <- dim_info[1]
  gene_count <- dim_info[2]

  if (sample_count <= gene_count) {
    stop(paste0(
      "\u26a0\ufe0f \u57fa\u56e0\u6570(", gene_count, ")\u4e0d\u5927\u4e8e\u6837\u672c\u6570(", sample_count,
      ")\uff0c\u505c\u6b62\u5206\u6790"
    ), call. = FALSE)
  }

  message(paste(
    "\u2705 \u57fa\u56e0\u6570(", gene_count, ") <\u6837\u672c\u6570(", sample_count,
    ")\uff0c\u7ee7\u7eed\u5206\u6790"
  ))

  # 存储表达矩阵
  seurat_obj@misc$MEscnet$datExpr <- merged_data

  # 初始化 TOM 列表
  TOMs <- list()

  tryCatch({
    # 运行 Mescnet 软阈值选择
    seurat_obj <- pick_power_scale_free_power(seurat_obj)

    # 设置 WGCNA 参数
    wgcna_name <- "MEscnet"
    seurat_obj <- hdWGCNA::SetActiveWGCNA(seurat_obj, wgcna_name)

    # 测试软阈值
    seurat_obj <- hdWGCNA::TestSoftPowers(
      seurat_obj,
      networkType = "signed"
    )

    # 构建网络
    seurat_obj <- hdWGCNA::ConstructNetwork(
      seurat_obj,
      tom_name = "hdwgcna"
    )

    # Mescnet 网络：pick_power_scale_free_power 中基于 pcor 的邻接矩阵
    TOM_mescnet <- seurat_obj@misc$MEscnet$tom
    if (is.null(TOM_mescnet)) {
      stop("Mescnet TOM (`seurat_obj@misc$MEscnet$tom`) is NULL; no soft-thresholding power was selected for the pcor network.")
    }
    TOMs[["MEscnet"]] <- TOM_mescnet

    # 获取 hdWGCNA 的 TOM 矩阵，并把基因顺序对齐到 Mescnet 网络
    TOM <- hdWGCNA::GetTOM(seurat_obj)
    colnames(TOM) <- rownames(TOM) <- colnames(TOM_mescnet)
    TOMs[["hdwgcna"]] <- TOM

    # 生成随机矩阵
    n <- length(colnames(TOM))
    random_adj_matrix <- matrix(0, nrow = n, ncol = n)

    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        random_adj_matrix[i, j] <- stats::runif(1, 0, 1)
        random_adj_matrix[j, i] <- random_adj_matrix[i, j]
      }
    }

    diag(random_adj_matrix) <- 0
    rownames(random_adj_matrix) <- colnames(random_adj_matrix) <- rownames(TOM)
    TOMs[["random"]] <- random_adj_matrix

    # 运行 benchmark 并绘图
    benchmark <- plot_metanet_auc(result = TOMs, cell_type = cell_type, go = go)

    if (save_TOMs) {
      save(TOMs, file = tom_file)
    }

    list(
      seurat_obj = seurat_obj,
      TOMs = TOMs,
      plot = benchmark$plot,
      summary_stats = benchmark$summary_stats,
      plot_df = benchmark$plot_df
    )

  }, error = function(e) {
    stop(paste("\u5206\u6790\u8fc7\u7a0b\u4e2d\u51fa\u9519:", e$message))
  })
}
