#' Aggregate spatial spots into metaspots (Seurat v5 compatible)
#'
#' Merges spatially adjacent spots within each level of `group.by` into
#' metaspots, the spatial analogue of metacells, to reduce sparsity in spatial
#' transcriptomics data.
#'
#' @param seurat_obj A Seurat object with spatial coordinates stored in
#'   `seurat_obj@meta.data` as `row`, `col`, `imagerow` and `imagecol`.
#' @param group.by Character vector of metadata columns defining the groups to
#'   aggregate within. Default `"seurat_clusters"`.
#' @param ident.group Metadata column used to set the identity of the resulting
#'   metaspot object; must be one of `group.by`. Default `"seurat_clusters"`.
#' @param assay Spatial assay to aggregate. Default `"Spatial"`.
#' @param slot,layer Layer or slot holding the expression values. `layer` is used
#'   with Seurat v5 and `slot` otherwise; both default to `"counts"`.
#' @param mode Aggregation mode, either `"sum"` (default) or `"average"`.
#' @param min_spots Minimum number of spots a group must have to be aggregated.
#'   Default `50`.
#' @param wgcna_name Name of the hdWGCNA experiment the metaspot object is
#'   stored under. `NULL` uses `seurat_obj@misc$active_wgcna`.
#'
#' @return The input Seurat object with the merged metaspot object stored under
#'   `wgcna_name`.
#'
#' @details
#' This function is a drop-in replacement for [hdWGCNA::MetaspotsByGroups()].
#' Compared with the hdWGCNA version it additionally:
#'
#' * writes the `ident.group` identity of each source group back into the
#'   metadata of the corresponding metaspot object, which is required for
#'   `Idents<-` to succeed on Seurat v5 objects;
#' * keeps `groupings` and the metaspot list in sync when a group fails to
#'   aggregate, instead of silently dropping only one of the two;
#' * joins assay layers after merging on Seurat v5 objects.
#'
#' @seealso [`MetacellsByGroups()`]
#'
#' @examples
#' \dontrun{
#' seurat_obj <- MetaspotsByGroups(
#'   seurat_obj,
#'   group.by = "cell_type",
#'   ident.group = "cell_type",
#'   assay = "Spatial",
#'   min_spots = 50
#' )
#' }
#'
#' @export
MetaspotsByGroups <- function(seurat_obj, group.by = c("seurat_clusters"), ident.group = "seurat_clusters",
                              assay = "Spatial", slot = "counts", layer = "counts", mode = "sum",
                              min_spots = 50, wgcna_name = NULL)
{
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  if (any(grepl("#", group.by))) {
    stop("Invalid character # found in group.by, please re-name the group.")
  }
  if (!(ident.group %in% group.by)) {
    stop("ident.group must be in group.by")
  }
  if (!(mode %in% c("sum", "average"))) {
    stop("Invalid choice for mode. Mode can be either sum or average.")
  }
  if (!all(c("row", "col", "imagerow", "imagecol") %in% colnames(seurat_obj@meta.data))) {
    stop("Spatial coordinates missing from seurat_obj@meta.data, must have columns named row, col, imagerow, and imagecol.")
  }
  if (is.null(assay)) {
    assay <- SeuratObject::DefaultAssay(seurat_obj)
  } else if (!(assay %in% names(seurat_obj@assays))) {
    stop(paste0("Assay ", assay, " not found in seurat_obj."))
  }

  if (length(group.by) > 1) {
    seurat_meta <- seurat_obj@meta.data[, group.by, drop = FALSE]
    for (col in colnames(seurat_meta)) {
      seurat_meta[[col]] <- as.character(seurat_meta[[col]])
    }
    seurat_obj$metacell_grouping <- apply(seurat_meta, 1, paste, collapse = "#")
  } else {
    seurat_obj$metacell_grouping <- as.character(seurat_obj@meta.data[[group.by]])
  }

  groupings <- unique(seurat_obj$metacell_grouping)
  groupings <- groupings[order(groupings)]
  group_counts <- table(seurat_obj$metacell_grouping) < min_spots

  if (any(group_counts)) {
    warning(paste0("Removing the following groups that did not meet min_spots: ",
                   paste(names(group_counts)[group_counts], collapse = ", ")))
  }
  groupings <- groupings[table(seurat_obj$metacell_grouping) >= min_spots]

  if (length(groupings) == 0) {
    stop("No groups met the min_spots requirement.")
  }

  seurat_list <- lapply(groupings, function(x) {
    seurat_obj[, seurat_obj$metacell_grouping == x]
  })
  names(seurat_list) <- groupings

  # 屏蔽底层烦人的警告
  metaspot_list <- suppressWarnings(mapply(hdWGCNA:::ConstructMetaspots, cur_seurat = seurat_list,
                                           MoreArgs = list(mode = mode, assay = assay, slot = slot, layer = layer)))
  names(metaspot_list) <- groupings

  remove <- which(sapply(metaspot_list, is.null))
  if (length(remove) > 0) {
    metaspot_list <- metaspot_list[-remove]
    groupings <- groupings[-remove] # 必须同步删减 groupings，防止错位
  }
  if (length(metaspot_list) == 0) {
    stop("All metaspot aggregations failed.")
  }

  # 🌟🌟🌟 核心绝杀修复区：手动将身份信息强行注入 Metadata 🌟🌟🌟
  for (k in seq_along(metaspot_list)) {
    # 从原始切片里找出这个组的一个代表细胞，提取它的真实 ident.group 值
    rep_cell <- colnames(seurat_list[[names(metaspot_list)[k]]])[1]
    ident_val <- as.character(seurat_obj@meta.data[rep_cell, ident.group])

    # 强行塞进新生成的 metaspot 对象里！
    metaspot_list[[k]]@meta.data[[ident.group]] <- ident_val
  }

  if (length(metaspot_list) > 1) {
    density_agg <- unlist(lapply(metaspot_list, function(x) x@misc$density_agg))
    density_orig <- unlist(lapply(metaspot_list, function(x) x@misc$density_orig))
    density_names <- names(density_agg)
  } else {
    density_agg <- metaspot_list[[1]]@misc$density_agg
    density_orig <- metaspot_list[[1]]@misc$density_orig
    density_names <- "all"
  }
  density_df <- data.frame(group = density_names, agg = as.numeric(density_agg), orig = as.numeric(density_orig))

  spot_neighbors <- list()
  for (k in 1:length(metaspot_list)) {
    spot_neighbors <- c(spot_neighbors, metaspot_list[[k]]@misc$spot_neighbors)
  }

  # 屏蔽 merge 产生的冗余警告，处理 V5 图层
  if (length(metaspot_list) > 1) {
    metaspot_obj <- suppressWarnings(merge(metaspot_list[[1]], metaspot_list[2:length(metaspot_list)]))
    if (as.numeric(substr(utils::packageVersion("Seurat"), 1, 1)) >= 5) {
      metaspot_obj[[assay]] <- suppressWarnings(SeuratObject::JoinLayers(metaspot_obj[[assay]]))
    }
  } else {
    metaspot_obj <- metaspot_list[[1]]
  }

  metaspot_obj@misc$spot_neighbors <- spot_neighbors
  metaspot_obj@misc$density <- density_df

  # 这里绝对不会再报错了！
  SeuratObject::Idents(metaspot_obj) <- metaspot_obj@meta.data[[ident.group]]
  seurat_obj <- hdWGCNA::SetMetacellObject(seurat_obj, metaspot_obj, wgcna_name)

  seurat_obj
}
