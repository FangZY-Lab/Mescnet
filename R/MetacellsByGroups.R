#' Aggregate single cells into metacells (Seurat v5 compatible)
#'
#' Transcriptionally neighbouring cells within each level of `group.by` are
#' aggregated into metacells, which reduces expression sparsity and stabilises
#' gene-gene association estimates while retaining local cellular heterogeneity.
#'
#' @param seurat_obj A Seurat object.
#' @param group.by Character vector of metadata columns defining the metacell
#'   groups. Default `"seurat_clusters"`.
#' @param ident.group Metadata column used to set the identity of the resulting
#'   metacell object; must be one of `group.by`. Default `"seurat_clusters"`.
#' @param k Number of nearest neighbours used to construct each metacell.
#'   Default `25`.
#' @param reduction Name of the dimensional reduction used to find neighbouring
#'   cells. Default `"pca"`.
#' @param dims Dimensions of `reduction` to use. `NULL` uses all dimensions
#'   stored in the reduction.
#' @param assay Assay to aggregate. `NULL` uses the default assay.
#' @param slot,data Layer or slot holding the expression values. `layer` is used
#'   with SeuratObject >= 5.0.0 and `slot` otherwise; both default to
#'   `"counts"`.
#' @param mode Aggregation mode, either `"average"` (default) or `"sum"`.
#' @param cells.use Optional vector of cells to restrict the aggregation to.
#' @param min_cells Minimum number of cells a group must have to be aggregated.
#'   Default `100`.
#' @param max_shared Maximum number of cells two metacells may share. Default
#'   `15`.
#' @param target_metacells Target number of metacells per group. Default
#'   `1000`.
#' @param max_iter Maximum number of iterations for the metacell construction.
#'   Default `5000`.
#' @param verbose Print progress messages. Default `FALSE`.
#' @param wgcna_name Name of the hdWGCNA experiment the metacell object is
#'   stored under. `NULL` uses `seurat_obj@misc$active_wgcna`.
#'
#' @return The input Seurat object with the merged metacell object and the
#'   metacell parameters stored under `wgcna_name`.
#'
#' @details
#' This function is a drop-in replacement for [hdWGCNA::MetacellsByGroups()] with
#' two differences:
#'
#' * the expression matrix is read with the `layer` argument when
#'   `SeuratObject >= 5.0.0` is installed and with `slot` otherwise, so it works
#'   on both Seurat v4 and v5 objects;
#' * after merging, layers are joined with [SeuratObject::JoinLayers()] on
#'   Seurat v5 objects, which is required before the metacell matrix can be
#'   normalised and used for network inference.
#'
#' Because it shares its name with the hdWGCNA function, load `Mescnet` after
#' `hdWGCNA` (or call `Mescnet::MetacellsByGroups()`) to make sure this version
#' is the one that gets used.
#'
#' @seealso [`CreaterMEscnet()`], [`MetaspotsByGroups()`]
#'
#' @examples
#' \dontrun{
#' seurat_obj <- MetacellsByGroups(
#'   seurat_obj,
#'   group.by = "cell_type",
#'   reduction = "harmony",
#'   k = 10,
#'   max_shared = 5,
#'   ident.group = "cell_type",
#'   target_metacells = 2000
#' )
#' }
#'
#' @export
MetacellsByGroups <- function(seurat_obj, group.by = c("seurat_clusters"), ident.group = "seurat_clusters",
    k = 25, reduction = "pca", dims = NULL, assay = NULL, slot = "counts",
    layer = "counts", mode = "average", cells.use = NULL, min_cells = 100,
    max_shared = 15, target_metacells = 1000, max_iter = 5000,
    verbose = FALSE, wgcna_name = NULL)
{
    if (is.null(wgcna_name)) {
        wgcna_name <- seurat_obj@misc$active_wgcna
    }
    hdWGCNA::CheckWGCNAName(seurat_obj, wgcna_name)
    if (any(grepl("#", group.by))) {
        stop("Invalid character # found in group.by, please re-name the group.")
    }
    if (!(ident.group %in% group.by)) {
        stop("ident.group must be in group.by")
    }
    if (!(mode %in% c("sum", "average"))) {
        stop("Invalid choice for mode. Mode can be either sum or average.")
    }
    if (!(reduction %in% names(seurat_obj@reductions))) {
        stop(paste0("Invalid reduction (", reduction, "). Reductions in Seurat object: ",
            paste(names(seurat_obj@reductions), collapse = ", ")))
    }
    if (is.null(assay)) {
        assay <- SeuratObject::DefaultAssay(seurat_obj)
    } else if (!(assay %in% names(seurat_obj@assays))) {
        stop(paste0("Assay ", assay, " not found in seurat_obj. Select a valid assay: ",
            paste0(names(seurat_obj@assays), collapse = ", ")))
    }

    # ============== 【核心修改区：让代码懂 Seurat V5 的 layer】 ==============
    if (!(layer %in% c("counts", "data", "scale.data")) && !(slot %in% c("counts", "data", "scale.data"))) {
        stop("Invalid input for layer/slot. Valid choices are counts, data, scale.data.")
    } else {
        # 动态判断：如果是 V5 就乖乖用 layer，否则用 slot
        if (utils::packageVersion("SeuratObject") >= "5.0.0") {
            slot_dim <- dim(SeuratObject::GetAssayData(seurat_obj, assay = assay, layer = layer))
        } else {
            slot_dim <- dim(SeuratObject::GetAssayData(seurat_obj, assay = assay, slot = slot))
        }

        if (is.null(slot_dim) || any(slot_dim == 0)) {
            stop(paste(c("Selected layer/slot not found in this assay.")))
        }
    }
    # =========================================================================

    if (min_cells < k) {
        warning("min_cells is smaller than k, this may result in downstream errors if very small groups are allowed.")
    }
    if (max_shared < 0) {
        warning(paste0("max_shared specified (", max_shared,
            ") is too low, setting max_shared <- 0"))
        max_shared <- 0
    }
    if (!is.null(cells.use)) {
        seurat_full <- seurat_obj
        seurat_obj <- seurat_obj[, cells.use]
    }
    if (length(group.by) > 1) {
        seurat_meta <- seurat_obj@meta.data[, group.by]
        for (col in colnames(seurat_meta)) {
            seurat_meta[[col]] <- as.character(seurat_meta[[col]])
        }
        seurat_obj$metacell_grouping <- apply(seurat_meta, 1,
            paste, collapse = "#")
    }
    else {
        seurat_obj$metacell_grouping <- as.character(seurat_obj@meta.data[[group.by]])
    }
    groupings <- unique(seurat_obj$metacell_grouping)
    groupings <- groupings[order(groupings)]
    group_counts <- table(seurat_obj$metacell_grouping) < min_cells
    if (any(group_counts)) {
        warning(paste0("Removing the following groups that did not meet min_cells: ",
            paste(names(group_counts)[group_counts], collapse = ", ")))
    }
    groupings <- groupings[table(seurat_obj$metacell_grouping) >= min_cells]
    if (length(groupings) == 0) {
        stop("No groups met the min_cells requirement.")
    }
    meta_df <- as.data.frame(do.call(rbind, strsplit(groupings, "#")))
    colnames(meta_df) <- group.by
    meta_list <- lapply(1:nrow(meta_df), function(i) {
        x <- list(as.character(meta_df[i, ]))[[1]]
        names(x) <- colnames(meta_df)
        x
    })
    seurat_list <- lapply(groupings, function(x) {
        seurat_obj[, seurat_obj$metacell_grouping == x]
    })
    names(seurat_list) <- groupings

    metacell_list <- mapply(hdWGCNA:::ConstructMetacells, seurat_obj = seurat_list,
        name = groupings, meta = meta_list, MoreArgs = list(k = k,
            reduction = reduction, dims = dims, assay = assay,
            slot = slot, layer = layer, return_metacell = TRUE,
            mode = mode, max_shared = max_shared, max_iter = max_iter,
            target_metacells = target_metacells, verbose = verbose,
            wgcna_name = wgcna_name))
    names(metacell_list) <- groupings
    remove <- which(sapply(metacell_list, is.null))
    if (length(remove) > 0) {
        metacell_list <- metacell_list[-remove]
    }
    if (length(metacell_list) == 0) {
        stop("Metacell construction failed for every group.")
    }
    run_stats <- as.data.frame(do.call(rbind, lapply(metacell_list,
        function(x) {
            x@misc$run_stats
        })))
    rownames(run_stats) <- 1:nrow(run_stats)
    for (i in 1:length(group.by)) {
        run_stats[[group.by[i]]] <- do.call(rbind, strsplit(as.character(run_stats$name),
            "#"))[, i]
    }
    if (length(metacell_list) > 1) {
        metacell_obj <- merge(metacell_list[[1]], metacell_list[2:length(metacell_list)])
        if (hdWGCNA::CheckSeurat5()) {
            metacell_obj <- SeuratObject::JoinLayers(metacell_obj)
        }
    }
    else {
        metacell_obj <- metacell_list[[1]]
    }
    SeuratObject::Idents(metacell_obj) <- metacell_obj@meta.data[[ident.group]]
    if (!is.null(cells.use)) {
        seurat_obj <- seurat_full
    }
    seurat_obj <- hdWGCNA::SetMetacellObject(seurat_obj, metacell_obj,
        wgcna_name)
    param_list <- list(metacell_k = k, metacell_reduction = reduction,
        metacell_assay = assay, metacell_stats = run_stats)
    if (hdWGCNA::CheckSeurat5()) {
        param_list["metacell_layer"] <- layer
    }
    else {
        param_list["metacell_slot"] <- slot
    }
    seurat_obj <- hdWGCNA::SetWGCNAParams(seurat_obj, params = param_list,
        wgcna_name)
    seurat_obj
}

