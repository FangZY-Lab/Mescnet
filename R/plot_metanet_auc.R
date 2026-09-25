#' Benchmark gene networks by EGAD neighbour-voting AUROC
#'
#' Scores each network in `result` by running Gene Ontology guided gene set
#' analysis (GBA) with [EGAD::run_GBA()] and plots the resulting neighbour
#' voting AUROC distributions as violin + notched boxplots.
#'
#' @param result A named list of square gene-gene association (or topological
#'   overlap) matrices, with gene symbols on the row and column names. Names are
#'   used as the network labels on the x axis, for example
#'   `list(MEscnet = ..., hdwgcna = ..., random = ...)`.
#' @param cell_type Label used for the title of the saved figure and for the
#'   output file name. Default `"test"`.
#' @param go Data frame of GO annotations with the GO term in column 2 and the
#'   gene symbol in column 3. Defaults to a `go` object found in the global
#'   environment (for backwards compatibility with the original scripts).
#' @param out_dir Directory the figure is written to. Pass `NULL` to skip
#'   writing a file. Default `"."`.
#' @param width,height Figure size in inches passed to [ggplot2::ggsave()].
#'   Defaults `8` and `6`.
#'
#' @return A list with
#'
#' * `plot` - the `ggplot` object.
#' * `summary_stats` - per-network number of GO terms, minimum, maximum and
#'   standard deviation of the AUROC values.
#' * `plot_df` - long data frame of every AUROC value with its network label.
#'
#' As a side effect the per-network summary statistics are printed.
#'
#' @details
#' Genes for which no Entrez ID can be mapped are dropped before the GBA run,
#' and the lower triangle of each matrix is restored from the upper triangle in
#' case it was zeroed out upstream.
#'
#' @seealso [`run_wgcna_analysis()`]
#'
#' @examples
#' \dontrun{
#' load("TOMs.Rdata")
#' res <- plot_metanet_auc(TOMs, cell_type = "INH", go = go)
#' res$summary_stats
#' }
#'
#' @export
plot_metanet_auc <- function(result, cell_type = "test", go = .Mescnet_default_go(),
                             out_dir = ".", width = 8, height = 6) {

  .Mescnet_require("AnnotationDbi", "org.Hs.eg.db", "EGAD")

  # 初始化结果数据框
  plot_df <- data.frame()

  # 遍历 result 中的每个矩阵
  for (i in 1:length(result)) {
    tom <- result[[i]]
    gene1 <- rownames(tom)

    # 将基因符号转换为 Entrez IDs
    gene_symbols <- AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = gene1,
      column = "ENTREZID",
      keytype = "SYMBOL",
      multiVals = "first"
    )

    # 创建基因 ID 和符号的数据框
    result_df <- data.frame(
      EntrezID = gene_symbols[gene1],
      Symbol = gene1
    )
    gene_table <- stats::na.omit(result_df)

    # 检查矩阵是否对称
    is_symmetric <- all(tom[upper.tri(tom)] == t(tom)[upper.tri(tom)])
    cat("Matrix", names(result)[i], "is symmetric:", is_symmetric, "\n")

    # 如果下三角被清零，用上三角恢复
    tom[lower.tri(tom)] <- t(tom)[lower.tri(tom)]

    # 筛选基因
    genelist <- rownames(tom)
    genelist <- genelist[(genelist %in% gene_table$Symbol)]
    tom <- tom[genelist, genelist]

    # 获取基因 ID
    ix <- match(genelist, gene_table$Symbol)
    gene_ids <- gene_table$EntrezID[ix]

    # 重命名 TOM 矩阵的行列
    colnames(tom) <- gene_ids
    rownames(tom) <- gene_ids

    goterms <- unique(go[, 3])
    annotations <- EGAD::make_annotations(go[, c(2, 3)], gene_ids, goterms)

    # 运行 GBA 分析
    GO_groups_voted <- EGAD::run_GBA(tom, annotations, max = Inf)

    # 提取 AUROC 值
    auc_GO_nv <- GO_groups_voted[[1]][, 1]
    auc_GO_nd <- GO_groups_voted[[1]][, 3]

    # 创建数据框
    df <- data.frame(
      nv_auc = auc_GO_nv,
      nd_auc = auc_GO_nd,
      tom = names(result)[i]
    )

    # 合并到 plot_df
    plot_df <- rbind(plot_df, df)
  }

  # 数据分组统计
  summary_stats <- plot_df %>%
    dplyr::group_by(tom) %>%
    dplyr::summarise(
      n = dplyr::n(),
      min_auc = min(nv_auc),
      max_auc = max(nv_auc),
      sd_auc = stats::sd(nv_auc),
      .groups = "drop"
    )
  print(summary_stats)

  # 确保 nv_auc 是数值型
  plot_df$nv_auc <- as.numeric(plot_df$nv_auc)

  # 绘制小提琴图和箱线图
  p1 <- ggplot2::ggplot(plot_df, ggplot2::aes(x = reorder(tom, nv_auc, median), y = nv_auc, fill = tom)) +
    ggplot2::geom_violin(trim = FALSE, scale = "width") +
    ggplot2::geom_boxplot(notch = TRUE, width = 0.25, fill = "white", outlier.size = 0) +
    ggplot2::ylim(0, 1) +
    ggplot2::xlab("") +
    ggplot2::ylab("AUC") +
    viridis::scale_fill_viridis(discrete = TRUE, option = "viridis") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none",
      axis.line = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title = ggplot2::element_text(vjust = -0.5, hjust = 0.5, face = "plain"),
      panel.grid.major.y = ggplot2::element_line(color = "lightgrey", linewidth = 0.5)
    ) +
    ggplot2::ggtitle("Neighbor voting")

  # 保存图
  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    ggplot2::ggsave(file.path(out_dir, paste0(cell_type, ".png")), p1, width = width, height = height)
  }

  # 返回结果数据框和统计信息
  list(plot = p1, summary_stats = summary_stats, plot_df = plot_df)
}


# 默认的 GO 注释表：回退到全局环境中的 `go` 对象，以兼容旧脚本。
.Mescnet_default_go <- function() {
  if (exists("go", envir = globalenv(), inherits = FALSE)) {
    return(get("go", envir = globalenv(), inherits = FALSE))
  }
  stop(
    "No `go` annotation table supplied and no `go` object found in the global ",
    "environment. Pass it explicitly, e.g. `plot_metanet_auc(TOMs, go = go)`.",
    call. = FALSE
  )
}
