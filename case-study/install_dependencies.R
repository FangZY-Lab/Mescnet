# Install the dependencies used by the PBMC case study.

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat")
}

remotes::install_github(
  "smorabit/hdWGCNA",
  dependencies = TRUE,
  upgrade = "never"
)

# Optional: only needed for Plot_MEscnet_Network().
remotes::install_github(
  "Jiawang1209/ggNetView",
  dependencies = TRUE,
  upgrade = "never"
)

# Optional: only needed for Plot_Module_Scores().
BiocManager::install("UCell", ask = FALSE, update = FALSE)
remotes::install_github(
  "FangZY-Lab/Mescnet",
  dependencies = TRUE,
  upgrade = "never"
)

message("Dependencies installed. Restart R before running the case study.")
