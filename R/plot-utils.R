# Shared palette and plotting helpers.
.Mescnet_palette <- c(
  "#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF",
  "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF",
  "#7E6148FF", "#B09C85FF", "#3B4992FF", "#EE0000FF",
  "#008B45FF", "#631879FF", "#008280FF", "#BB0021FF",
  "#5F559BFF", "#A20056FF", "#808180FF", "#1B1919FF"
)

#' Create a consistent Mescnet colour palette
#'
#' Returns a fixed high-contrast palette for cell groups and module colours.
#' When more colours are requested than the palette contains, the colours are
#' interpolated.
#'
#' @param n Number of colours to return.
#'
#' @return A character vector of hexadecimal colour values.
#' @export
#'
#' @examples
#' Mescnet_colors(4)
Mescnet_colors <- function(n = length(.Mescnet_palette)) {
  if (length(n) != 1L || is.na(n) || n < 0 || n != as.integer(n)) {
    stop("`n` must be a single non-negative integer.", call. = FALSE)
  }
  n <- as.integer(n)
  if (n == 0L) {
    return(character())
  }
  if (n <= length(.Mescnet_palette)) {
    return(.Mescnet_palette[seq_len(n)])
  }
  grDevices::colorRampPalette(.Mescnet_palette)(n)
}

.Mescnet_assert_seurat <- function(x, arg = "seurat_obj") {
  if (!inherits(x, "Seurat")) {
    stop("`", arg, "` must be a Seurat object.", call. = FALSE)
  }
  invisible(x)
}

.Mescnet_get_layer <- function(object, assay, layer = "data") {
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    return(Seurat::GetAssayData(object, assay = assay, layer = layer))
  }
  Seurat::GetAssayData(object, assay = assay, slot = layer)
}

.Mescnet_output_dir <- function(path) {
  if (is.null(path)) {
    return(invisible(NULL))
  }
  directory <- dirname(path)
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(NULL)
}

.Mescnet_module_names <- function(module_genes) {
  if (!is.list(module_genes) || !length(module_genes)) {
    stop("`module_genes` must be a non-empty list of gene vectors.",
         call. = FALSE)
  }
  module_names <- names(module_genes)
  if (is.null(module_names) || any(!nzchar(module_names))) {
    module_names <- paste0("Module", seq_along(module_genes))
  }
  names(module_genes) <- make.unique(module_names)
  module_genes
}

.Mescnet_complete_palette <- function(palette, levels) {
  levels <- as.character(levels)
  if (!length(levels)) {
    return(character())
  }
  if (is.null(palette) || !length(palette)) {
    colors <- Mescnet_colors(length(levels))
    names(colors) <- levels
    return(colors)
  }
  if (is.null(names(palette)) || any(!nzchar(names(palette)))) {
    colors <- rep_len(as.character(palette), length(levels))
    names(colors) <- levels
    return(colors)
  }
  colors <- as.character(palette[levels])
  missing <- is.na(colors) | !nzchar(colors)
  if (any(missing)) {
    colors[missing] <- Mescnet_colors(sum(missing))
  }
  names(colors) <- levels
  colors
}

.Mescnet_check_top_n <- function(top_n) {
  if (length(top_n) != 1L || is.na(top_n) || top_n < 1 ||
      top_n != as.integer(top_n)) {
    stop("`top_n` must be a single positive integer.", call. = FALSE)
  }
  as.integer(top_n)
}

.Mescnet_module_kme <- function(seurat_obj, tom_matrix, module_genes,
                                assay = "RNA") {
  module_genes <- .Mescnet_module_names(module_genes)
  tom_matrix <- as.matrix(tom_matrix)
  if (is.null(rownames(tom_matrix))) {
    stop("`tom_matrix` must have gene names in its row names.", call. = FALSE)
  }

  expr_data <- as.matrix(.Mescnet_get_layer(seurat_obj, assay, "data"))
  genes <- intersect(rownames(tom_matrix), rownames(expr_data))
  if (!length(genes)) {
    stop("No genes are shared by `tom_matrix` and the assay expression data.",
         call. = FALSE)
  }
  expr_data <- expr_data[genes, , drop = FALSE]
  module_names <- names(module_genes)

  module_eigengenes <- matrix(
    0,
    nrow = ncol(expr_data),
    ncol = length(module_genes),
    dimnames = list(colnames(expr_data), module_names)
  )
  for (i in seq_along(module_genes)) {
    valid_genes <- intersect(module_genes[[i]], genes)
    if (length(valid_genes)) {
      module_eigengenes[, i] <- colMeans(
        expr_data[valid_genes, , drop = FALSE]
      )
    }
  }

  kme <- as.data.frame(stats::cor(
    t(expr_data),
    module_eigengenes,
    use = "pairwise.complete.obs"
  ))
  colnames(kme) <- module_names

  list(
    expr_data = expr_data,
    genes = genes,
    kme = kme,
    module_genes = module_genes,
    module_names = module_names
  )
}

.Mescnet_hub_genes <- function(kme, module_genes, genes, top_n) {
  top_n <- .Mescnet_check_top_n(top_n)
  module_genes <- .Mescnet_module_names(module_genes)
  hub_genes <- character()
  for (i in seq_along(module_genes)) {
    module_name <- names(module_genes)[i]
    module_genes_i <- intersect(module_genes[[i]], genes)
    if (!length(module_genes_i)) {
      next
    }
    values <- abs(kme[module_genes_i, module_name])
    names(values) <- module_genes_i
    values <- values[is.finite(values)]
    if (!length(values)) {
      next
    }
    hub_genes <- c(
      hub_genes,
      names(utils::head(sort(values, decreasing = TRUE), top_n))
    )
  }
  unique(hub_genes)
}
