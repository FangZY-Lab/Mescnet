#' Create a Mescnet network object
#'
#' `CreaterMEscnet()` is the main entry point of the package. It selects the
#' genes used for network reconstruction, aggregates transcriptionally
#' neighbouring cells into metacells within each level of `group.by`, normalises
#' the metacell expression matrix, and stores it as `DatExpr` for the cell group
#' of interest.
#'
#' All results are written into the Seurat object rather than returned
#' separately:
#'
#' * `seurat_obj@misc$MEscnet$datExpr` - metacell expression matrix used for
#'   network inference.
#' * `seurat_obj@misc$MEscnet$MEscnet_genes` - genes selected for the network.
#' * `seurat_obj@misc$MEscnet$MEscnet_params` - parameters used for metacell
#'   aggregation.
#' * `seurat_obj@misc$MEscnet$MEscnet_metacell_obj` - the metacell Seurat object.
#' * `seurat_obj@misc$active_MEscnet` - name of the active network
#'   (`"MEscnet"`).
#'
#' @param seurat_obj A Seurat object. Must contain the dimensionality reduction
#'   named in `reduction`.
#' @param nfeatures Number of highly variable genes to use when `gene_list` is
#'   `NULL`. Default `2000`.
#' @param group.by Character vector of metadata columns used to define the
#'   metacell groups. Default `"cell_type"`.
#' @param reduction Name of the dimensional reduction used to find
#'   transcriptionally neighbouring cells. Default `"harmony"`.
#' @param k Number of nearest neighbours used to build each metacell.
#'   Default `10`.
#' @param max_shared Maximum number of cells that two metacells may share.
#'   Default `5`.
#' @param target_metacells Target number of metacells per group. Default `2000`.
#' @param group_name Level(s) of `group.by` used to build `DatExpr`. This
#'   selects which cell population the network is inferred for. Default
#'   `"INH"`.
#' @param gene_list Optional character vector of genes to use instead of
#'   computing highly variable genes. Default `NULL`.
#'
#' @return The input Seurat object with metacells, `DatExpr` and the network
#'   metadata added to `@misc$MEscnet`.
#'
#' @details
#' `CreaterMEscnet()` normalises the metacell expression matrix before
#' `DatExpr` is set. If you need the raw metacell counts instead (for example
#' when the metacell matrix is later used with tools that expect counts), use
#' [`CreaterMEscnet_data()`].
#'
#' @seealso [`CreaterMEscnet_data()`], [`pick_power_scale_free_power()`],
#'   [`ComputeMEscnetModules()`]
#'
#' @examples
#' \dontrun{
#' seurat_obj <- CreaterMEscnet(
#'   seurat_obj,
#'   nfeatures = 1000,
#'   group.by = "All_Cells",
#'   reduction = "umap",
#'   k = 10,
#'   max_shared = 10,
#'   target_metacells = 500,
#'   group_name = "All_Immune",
#'   gene_list = common_genes
#' )
#'
#' seurat_obj <- pick_power_scale_free_power(seurat_obj)
#' seurat_obj <- ComputeMEscnetModules(seurat_obj, resolution = 1, number = 1)
#' }
#'
#' @export
CreaterMEscnet <- function(
    seurat_obj,
    nfeatures = 2000,
    group.by = "cell_type",
    reduction = "harmony",
    k = 10,
    max_shared = 5,
    target_metacells = 2000,
    group_name = "INH",
    gene_list = NULL
) {

  if (!is.null(gene_list)) {
    # 如果提供了gene_list，使用它设置WGCNA
    seurat_obj <- hdWGCNA::SetupForWGCNA(
      seurat_obj,
      gene_select = "custom",
      gene_list = gene_list,
      wgcna_name = "MEscnet"
    )
    message("WGCNA set up with provided gene_list.")
  } else {
    # 如果未提供gene_list，计算高变基因
    seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = nfeatures)
    var_genes <- Seurat::VariableFeatures(seurat_obj)
    seurat_obj <- hdWGCNA::SetupForWGCNA(
      seurat_obj,
      gene_select = "custom",
      gene_list = var_genes,
      wgcna_name = "MEscnet"
    )
    message("WGCNA set up with ", length(var_genes), " variable genes.")
  }

  # 通过 group.by 创建 metacells
  seurat_obj <- MetacellsByGroups(
    seurat_obj = seurat_obj,
    group.by = group.by,
    reduction = reduction,
    k = k,
    max_shared = max_shared,
    ident.group = group.by,
    target_metacells = target_metacells
  )

  # 归一化 metacell 表达矩阵
  seurat_obj <- hdWGCNA::NormalizeMetacells(seurat_obj)

  # 设置 DatExpr
  seurat_obj <- hdWGCNA::SetDatExpr(
    seurat_obj,
    group_name = group_name,
    group.by = group.by,
    assay = "RNA",
    layer = "data"
  )

  seurat_obj <- .Mescnet_rename_wgcna_slots(seurat_obj)

  seurat_obj
}


#' Create a Mescnet network object from metacell counts
#'
#' `CreaterMEscnet_data()` mirrors [`CreaterMEscnet()`] but keeps the raw
#' metacell count matrix instead of normalising it. Metacells are aggregated
#' from the `"data"` layer, the aggregated matrix is copied into the `"data"`
#' layer of the metacell object, and `DatExpr` is set from those counts.
#'
#' Use this variant when the downstream network should be inferred on
#' non-normalised metacell counts; use [`CreaterMEscnet()`] for the standard
#' normalised workflow.
#'
#' @inheritParams CreaterMEscnet
#'
#' @return The input Seurat object with metacells, `DatExpr` and the network
#'   metadata added to `@misc$MEscnet`.
#'
#' @seealso [`CreaterMEscnet()`]
#'
#' @examples
#' \dontrun{
#' seurat_obj <- CreaterMEscnet_data(seurat_obj, group_name = "INH")
#' }
#'
#' @export
CreaterMEscnet_data <- function(
    seurat_obj,
    nfeatures = 2000,
    group.by = "cell_type",
    reduction = "harmony",
    k = 10,
    max_shared = 5,
    target_metacells = 2000,
    group_name = "INH",
    gene_list = NULL
) {

  if (!is.null(gene_list)) {
    seurat_obj <- hdWGCNA::SetupForWGCNA(
      seurat_obj,
      gene_select = "custom",
      gene_list = gene_list,
      wgcna_name = "MEscnet"
    )
    message("WGCNA set up with provided gene_list.")
  } else {
    seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = nfeatures)
    var_genes <- Seurat::VariableFeatures(seurat_obj)
    seurat_obj <- hdWGCNA::SetupForWGCNA(
      seurat_obj,
      gene_select = "custom",
      gene_list = var_genes,
      wgcna_name = "MEscnet"
    )
    message("WGCNA set up with ", length(var_genes), " variable genes.")
  }

  # 通过 group.by 创建 metacells（使用 data 层）
  seurat_obj <- MetacellsByGroups(
    seurat_obj = seurat_obj,
    group.by = group.by,
    reduction = reduction,
    k = k,
    layer = "data",
    max_shared = max_shared,
    ident.group = group.by,
    target_metacells = target_metacells
  )

  # 用 counts 覆盖 data 层，使 DatExpr 建立在原始 counts 上
  seurat_obj@misc$MEscnet$wgcna_metacell_obj@assays$RNA$data <-
    seurat_obj@misc$MEscnet$wgcna_metacell_obj@assays$RNA$counts

  # 设置 DatExpr
  seurat_obj <- hdWGCNA::SetDatExpr(
    seurat_obj,
    group_name = group_name,
    group.by = group.by,
    assay = "RNA",
    layer = "data"
  )

  seurat_obj <- .Mescnet_rename_wgcna_slots(seurat_obj)

  seurat_obj
}


# ---------------------------------------------------------------------------
# Internal helper: move hdWGCNA's slot names over to the Mescnet namespace.
# ---------------------------------------------------------------------------
.Mescnet_rename_wgcna_slots <- function(seurat_obj) {

  # 确保 active network 指向 MEscnet
  if ("active_wgcna" %in% names(seurat_obj@misc)) {
    seurat_obj@misc$active_MEscnet <- seurat_obj@misc$active_wgcna
    seurat_obj@misc$active_wgcna <- NULL
  } else {
    message("active_wgcna not exist")
  }

  # 确保 MEscnet 存在
  if (!("MEscnet" %in% names(seurat_obj@misc))) {
    message("MEscnet 不存在于 Seurat 对象的 @misc 列表中")
    return(seurat_obj)
  }

  # 依次重命名子元素
  renames <- c(
    wgcna_genes         = "MEscnet_genes",
    wgcna_params        = "MEscnet_params",
    wgcna_metacell_obj  = "MEscnet_metacell_obj"
  )

  for (old in names(renames)) {
    if (old %in% names(seurat_obj@misc$MEscnet)) {
      new <- renames[[old]]
      seurat_obj@misc$MEscnet[[new]] <- seurat_obj@misc$MEscnet[[old]]
      seurat_obj@misc$MEscnet[[old]] <- NULL
    }
  }

  seurat_obj
}
