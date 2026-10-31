# ==============================================================================
# Mescnet PBMC case study (Seurat v5)
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Mescnet)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

# Named colours are matched to cell-group labels by the plotting functions.
GLOBAL_NATURE_COLORS <- c(
  "CD14+ Mono" = "#4DBBD5FF",
  "FCGR3A+ Mono" = "#8491B4FF",
  "Platelet" = "#7E6148FF",
  "DC" = "#DC0000FF",
  "B" = "#3C5488FF",
  "CD8 T" = "#F39B7FFF",
  "NK" = "#91D1C2FF",
  "Naive CD4 T" = "#E64B35FF",
  "Memory CD4 T" = "#00A087FF"
)

data_file <- "pbmc.Rdata"
out_dir <- file.path("case-study", "figures")

data_url <- paste0(
  "https://github.com/FangZY-Lab/Mescnet/releases/download/",
  "untagged-d0b110f2c48eeb2bfd42/pbmc.Rdata"
)
if (!file.exists(data_file)) {
  message("Downloading the PBMC example data to ", data_file, " ...")
  utils::download.file(data_url, data_file, mode = "wb")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
load(data_file)
if (!exists("seurat_obj")) {
  stop("`pbmc.Rdata` must contain an object named `seurat_obj`.",
       call. = FALSE)
}

# ------------------------------------------------------------------------------
# 1. Build metacells and select the genes used for network inference
# ------------------------------------------------------------------------------
seurat_obj <- CreaterMEscnet(
  seurat_obj,
  nfeatures = 203,
  group.by = "cell_type",
  reduction = "umap",
  k = 10,
  max_shared = 5,
  target_metacells = 100,
  group_name = c("CD14+ Mono", "CD8 T", "B")
)

# ------------------------------------------------------------------------------
# 2. Estimate association matrices and soft-thresholding powers
# ------------------------------------------------------------------------------
seurat_obj <- pick_power_scale_free_power(seurat_obj)

# ------------------------------------------------------------------------------
# 3. Detect modules. number = 3 selects the bayes_cor2 matrix.
# ------------------------------------------------------------------------------
seurat_obj <- ComputeMEscnetModules(
  seurat_obj,
  min_genes_per_module = 18,
  resolution = 0.5,
  number = 3
)

module_genes1 <- seurat_obj@misc$MEscnet$MEscnet_modules$gene_lists
tom_matrix <- seurat_obj@misc$MEscnet$bayes_cor2

# ------------------------------------------------------------------------------
# 4. Plot the case-study results
# ------------------------------------------------------------------------------
Plot_UMAP(
  seurat_obj,
  group.by = "cell_type",
  save_path = file.path(out_dir, "UMAP_PBMC_Nature.pdf"),
  palette = GLOBAL_NATURE_COLORS
)

Plot_Module_Scores(
  seurat_obj,
  module_genes = module_genes1,
  group.by = "cell_type",
  feature_file = file.path(out_dir, "FeaturePlot_Modules.pdf"),
  violin_file = file.path(out_dir, "VlnPlot_Modules.pdf"),
  palette = GLOBAL_NATURE_COLORS
)

Plot_MEscnet_Network(
  seurat_obj,
  tom_matrix = tom_matrix,
  module_genes = module_genes1,
  top_n = 5,
  output_file = file.path(out_dir, "MEscnet_ggNetView_Nature_Hubs.pdf")
)

Plot_Hub_Heatmap(
  seurat_obj,
  tom_matrix = tom_matrix,
  module_genes = module_genes1,
  target_cells = c(
    "CD14+ Mono", "FCGR3A+ Mono", "DC",
    "B", "NK", "CD8 T"
  ),
  group.by = "cell_type",
  top_n = 5,
  palette = GLOBAL_NATURE_COLORS,
  output_file = file.path(out_dir, "Figure_E_HubGenes_Heatmap.pdf")
)

message("Case-study figures written to: ", normalizePath(out_dir))
