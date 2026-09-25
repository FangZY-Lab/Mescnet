#' Select soft-thresholding powers for the Mescnet association matrices
#'
#' Builds seven candidate gene-gene association matrices from the metacell
#' expression matrix (`seurat_obj@misc$MEscnet$datExpr`), evaluates the
#' scale-free topology fit of each over `powerVector`, and records the first
#' power whose R-squared exceeds `RsquaredCut`.
#'
#' The seven matrices combine the two statistical views of gene-gene
#' relationships that motivate Mescnet: the marginal Pearson correlation
#' (`cor`) and the partial correlation (`pcor`), which estimates conditional
#' associations given all other genes.
#'
#' @param seurat_obj A Seurat object containing
#'   `seurat_obj@misc$MEscnet$datExpr`, as produced by [`CreaterMEscnet()`].
#' @param powerVector Integer vector of candidate soft-thresholding powers.
#'   Default `1:20`.
#' @param RsquaredCut R-squared threshold used to pick the recommended power.
#'   Default `0.80`.
#' @param nBreaks Number of breaks used to bin the connectivity distribution
#'   before fitting the log-log scale-free model. Default `10`.
#' @param plot_diagnostics Draw the four-panel diagnostic plot for each
#'   candidate matrix. Default `TRUE`; set to `FALSE` for headless or batch runs.
#'
#' @return The input Seurat object with the following added to
#'   `@misc$MEscnet`:
#'
#' * `MEscnet_ppcor` - the partial correlation matrix.
#' * `bayes_cor1` ... `bayes_cor6` - the candidate BayesCor association
#'   matrices.
#' * `result_tables` - per-matrix diagnostic table (R-squared, slope, mean,
#'   median and maximum connectivity for every power).
#' * `recommended_powers` - named list of the selected power per matrix, `NA`
#'   when no power reached `RsquaredCut`.
#' * `adj_matrices` - named list of soft-thresholded adjacency matrices for the
#'   matrices with a selected power.
#' * `tom` - the adjacency matrix of the `pcor` network, used by
#'   [`ComputeMEscnetModules()`] as the default network.
#'
#' When `plot_diagnostics = TRUE`, a four-panel diagnostic plot (R-squared, mean,
#' median and maximum connectivity versus power) is drawn for each candidate
#' matrix as a side effect.
#'
#' @details
#' The six BayesCor scores are defined as:
#'
#' * `bayes_cor1`: `abs(pcor) * abs(cor)`.
#' * `bayes_cor2`: `pmax(pcor, 0) * pmax(cor, 0)`.
#' * `bayes_cor3`: min-max normalised `pcor * cor`.
#' * `bayes_cor4`: signed-hybrid normalised `pcor * cor`, i.e. rescaling both
#'   matrices to `[0, 1]` with `(1 + x) / 2`.
#' * `bayes_cor5`: z-scored `pcor * cor`.
#' * `bayes_cor6`: `((1 + pmax(pcor, 0)) / 2) * pmax(cor, 0)`.
#'
#' `NA` values produced by zero-variance genes are replaced by `0`.
#'
#' @seealso [`CreaterMEscnet()`], [`ComputeMEscnetModules()`]
#'
#' @examples
#' \dontrun{
#' seurat_obj <- pick_power_scale_free_power(seurat_obj)
#' seurat_obj@misc$MEscnet$recommended_powers
#' }
#'
#' @export
pick_power_scale_free_power <- function(seurat_obj, powerVector = 1:20, RsquaredCut = 0.80,
                                        nBreaks = 10, plot_diagnostics = TRUE) {

  datExpr <- seurat_obj@misc$MEscnet$datExpr
  if (is.null(datExpr)) {
    stop("`seurat_obj@misc$MEscnet$datExpr` is NULL. Run `CreaterMEscnet()` first.", call. = FALSE)
  }

  # 计算偏相关矩阵
  pcor_matrix <- ppcor::pcor(datExpr)$estimate
  diag(pcor_matrix) <- 0
  colnames(pcor_matrix) <- colnames(datExpr)
  rownames(pcor_matrix) <- colnames(datExpr)
  seurat_obj@misc$MEscnet$MEscnet_ppcor <- as.matrix(pcor_matrix)

  # 计算 cor 矩阵
  cor_matrix <- stats::cor(datExpr, method = "pearson")
  diag(cor_matrix) <- 0
  colnames(cor_matrix) <- colnames(datExpr)
  rownames(cor_matrix) <- colnames(datExpr)

  # 计算 BayesCor1 到 BayesCor7
  # BayesCor1: abs(pcor) * abs(cor)
  bayes_cor1 <- abs(pcor_matrix) * abs(cor_matrix)
  bayes_cor1[is.na(bayes_cor1)] <- 0
  seurat_obj@misc$MEscnet$bayes_cor1 <- as.matrix(bayes_cor1)

  # BayesCor2: 小于 0 的值设为 0
  bayes_cor2 <- pmax(pcor_matrix, 0) * pmax(cor_matrix, 0)
  bayes_cor2[is.na(bayes_cor2)] <- 0
  seurat_obj@misc$MEscnet$bayes_cor2 <- as.matrix(bayes_cor2)

  # BayesCor3: Min-Max 归一化到 [0, 1]
  cor_norm <- (cor_matrix - min(cor_matrix, na.rm = TRUE)) /
    (max(cor_matrix, na.rm = TRUE) - min(cor_matrix, na.rm = TRUE))
  pcor_norm <- (pcor_matrix - min(pcor_matrix, na.rm = TRUE)) /
    (max(pcor_matrix, na.rm = TRUE) - min(pcor_matrix, na.rm = TRUE))
  bayes_cor3 <- pcor_norm * cor_norm
  bayes_cor3[is.na(bayes_cor3)] <- 0
  seurat_obj@misc$MEscnet$bayes_cor3 <- as.matrix(bayes_cor3)

  # BayesCor4: Signed Hybrid 归一化到 [0, 1]
  cor_norm <- (1 + cor_matrix) / 2
  pcor_norm <- (1 + pcor_matrix) / 2
  bayes_cor4 <- pcor_norm * cor_norm
  bayes_cor4[is.na(bayes_cor4)] <- 0
  seurat_obj@misc$MEscnet$bayes_cor4 <- as.matrix(bayes_cor4)

  # BayesCor5: Z 分数标准化
  cor_z <- (cor_matrix - mean(cor_matrix, na.rm = TRUE)) / stats::sd(cor_matrix, na.rm = TRUE)
  pcor_z <- (pcor_matrix - mean(pcor_matrix, na.rm = TRUE)) / stats::sd(pcor_matrix, na.rm = TRUE)
  bayes_cor5 <- pcor_z * cor_z
  bayes_cor5[is.na(bayes_cor5)] <- 0
  seurat_obj@misc$MEscnet$bayes_cor5 <- as.matrix(bayes_cor5)

  # BayesCor6: 新的混合方法
  cor_norm2 <- (1 + pmax(pcor_matrix, 0)) / 2
  bayes_cor6 <- cor_norm2 * pmax(cor_matrix, 0)
  bayes_cor6[is.na(bayes_cor6)] <- 0
  seurat_obj@misc$MEscnet$bayes_cor6 <- as.matrix(bayes_cor6)

  # 定义所有关联矩阵
  cor_matrices <- list(
    pcor = pcor_matrix,
    bayes_cor1 = bayes_cor1,
    bayes_cor2 = bayes_cor2,
    bayes_cor3 = bayes_cor3,
    bayes_cor4 = bayes_cor4,
    bayes_cor5 = bayes_cor5,
    bayes_cor6 = bayes_cor6
  )

  # 初始化结果列表
  result_tables <- list()
  recommended_powers <- list()
  adj_matrices <- list()

  # 为每种关联矩阵计算软阈值
  for (name in names(cor_matrices)) {
    result <- data.frame(
      Power = powerVector,
      Rsquared = NA,
      Slope = NA,
      MeanK = NA,
      MedianK = NA,
      MaxK = NA
    )

    # 遍历 power 值
    for (i in seq_along(powerVector)) {
      power <- powerVector[i]
      adj <- abs(cor_matrices[[name]])^power
      diag(adj) <- 0
      k <- rowSums(adj, na.rm = TRUE)

      # 计算连接度分布
      hist_k <- graphics::hist(k, breaks = nBreaks, plot = FALSE)
      counts <- hist_k$counts
      mids <- hist_k$mids
      valid <- (mids > 0) & (counts > 0)
      if (length(mids[valid]) > 1) {
        log_k <- log10(mids[valid])
        log_p <- log10(counts[valid] / sum(counts[valid]))
        fit <- stats::lm(log_p ~ log_k)
        result$Slope[i] <- stats::coef(fit)[2]
        result$Rsquared[i] <- summary(fit)$r.squared
      }
      result$MeanK[i] <- mean(k, na.rm = TRUE)
      result$MedianK[i] <- stats::median(k, na.rm = TRUE)
      result$MaxK[i] <- max(k, na.rm = TRUE)
    }

    # 存储结果
    result_tables[[name]] <- result

    # 选择推荐的 power
    idx_r2 <- which(result$Rsquared > RsquaredCut)[1]
    power_selected <- ifelse(is.na(idx_r2), NA, result$Power[idx_r2])
    recommended_powers[[name]] <- power_selected

    # 计算邻接矩阵
    if (!is.na(power_selected)) {
      adj <- abs(cor_matrices[[name]])^power_selected
      diag(adj) <- 0
      adj_matrices[[name]] <- adj
    }

    # 绘制诊断图
    if (plot_diagnostics) {
      graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
      colors <- c("#00a6ac", "#45b97c", "#009ad6", "#769149")

      # 图1：R² vs power
      graphics::plot(result$Power, result$Rsquared, type = "p", pch = 19, cex = 2,
           col = colors[1], xlab = "Soft Threshold (power)",
           ylab = "Scale Free Topology Fit (R\u00b2)", main = paste(name, "R\u00b2 vs Power"))
      graphics::abline(h = RsquaredCut, col = "red", lty = 2)
      if (!is.na(power_selected)) {
        graphics::points(power_selected, result$Rsquared[idx_r2], pch = 19, cex = 2.5, col = "red")
        graphics::text(power_selected, result$Rsquared[idx_r2],
             labels = paste0("Power=", power_selected), pos = 3, cex = 0.9, col = "red")
      }

      # 图2：MeanK vs power
      graphics::plot(result$Power, result$MeanK, type = "p", pch = 19, cex = 2,
           col = colors[2], xlab = "Soft Threshold (power)",
           ylab = "Mean Connectivity", main = paste(name, "Mean Connectivity"))

      # 图3：MedianK vs power
      graphics::plot(result$Power, result$MedianK, type = "p", pch = 19, cex = 2,
           col = colors[3], xlab = "Soft Threshold (power)",
           ylab = "Median Connectivity", main = paste(name, "Median Connectivity"))

      # 图4：MaxK vs power
      graphics::plot(result$Power, result$MaxK, type = "p", pch = 19, cex = 2,
           col = colors[4], xlab = "Soft Threshold (power)",
           ylab = "Max Connectivity", main = paste(name, "Max Connectivity"))
    }

    cat("\u2705 ", name, "Recommended power =", power_selected, "\n")
  }

  # 存储结果
  seurat_obj@misc$MEscnet$result_tables <- result_tables
  seurat_obj@misc$MEscnet$recommended_powers <- recommended_powers
  seurat_obj@misc$MEscnet$adj_matrices <- adj_matrices
  seurat_obj@misc$MEscnet$tom <- adj_matrices[["pcor"]]  # 保留原始 pcor 的 TOM

  seurat_obj
}
