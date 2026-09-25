# Mescnet

**Mescnet** (*maximum-entropy single-cell gene network*) reconstructs
context-specific gene programs from single-cell and spatial transcriptomics
data.

Cells are first aggregated into local **metacells**, which reduces expression
sparsity and stabilises gene-gene association estimates while retaining local
heterogeneity. Gene relationships are then estimated from two complementary
statistical views:

- **marginal correlation**, which captures coordinated variation between a gene
  pair, and
- **partial correlation**, which estimates the conditional association after
  accounting for all other genes and follows from a maximum-entropy
  formulation of the expression distribution.

The two signals are combined into a `BayesCor` score, soft-thresholded into a
weighted network, and gene communities are detected with the Leiden algorithm.
The resulting modules can be projected back onto cells and spatial coordinates
as phenotype-associated, context-specific, spatially organised or dynamically
rewired gene programs.

Analysis code for the manuscript figures is in
[Mescnet_paper_code](https://github.com/Fanglabyaominghui/Mescnet_paper_code).

## Installation

```r
# hdWGCNA and its dependencies
remotes::install_github("smorabit/hdWGCNA")

# Mescnet
remotes::install_github("Fanglabyaominghui/Mescnet")

# optional, only needed for the EGAD benchmark functions
install.packages("EGAD")
BiocManager::install(c("AnnotationDbi", "org.Hs.eg.db"))
```

`Mescnet` requires R >= 4.1, Seurat v5 (`SeuratObject` >= 5.0.0) and
`igraph` >= 1.3.0 for `cluster_leiden()`.

## Quick start

```r
library(Seurat)
library(hdWGCNA)
library(Mescnet)   # load after hdWGCNA so Mescnet's MetacellsByGroups wins

seurat_obj <- CreaterMEscnet(
  seurat_obj,
  group.by        = "cell_type",
  reduction       = "harmony",
  k               = 10,
  max_shared      = 5,
  target_metacells = 2000,
  group_name      = "INH",
  gene_list       = genes
)

# partial-correlation network, then the six BayesCor variants
seurat_obj <- pick_power_scale_free_power(seurat_obj)

# Leiden modules from the pcor network (1); 2-7 select bayes_cor1-6
seurat_obj <- ComputeMEscnetModules(seurat_obj, min_genes_per_module = 20,
                                    resolution = 1, number = 1)
```

## What ends up in the object

Everything is stored under `seurat_obj@misc$MEscnet`:

| Element | Content |
| --- | --- |
| `datExpr` | metacell expression matrix, cells x genes |
| `MEscnet_genes` | genes selected for the network |
| `MEscnet_params` | metacell aggregation parameters |
| `MEscnet_metacell_obj` | the metacell Seurat object |
| `MEscnet_ppcor` | partial correlation matrix |
| `bayes_cor1` ... `bayes_cor6` | integrated association matrices |
| `result_tables` | scale-free fit diagnostics per power |
| `recommended_powers` | selected power per association matrix |
| `adj_matrices` | soft-thresholded adjacency matrices |
| `tom` | adjacency matrix of the `pcor` network |
| `moduleColors` | module label per gene, `"grey"` for unassigned |
| `MEscnet_modules` | `module_ids`, `module_sizes`, `gene_lists` |
| `MEscnet_network` | the clustered `igraph` graph |

The active network is tracked in `seurat_obj@misc$active_MEscnet`.

```r
datExpr <- seurat_obj@misc$MEscnet$datExpr
colors  <- seurat_obj@misc$MEscnet$moduleColors
modules <- seurat_obj@misc$MEscnet$MEscnet_modules$gene_lists
```

## Choosing an association matrix

`ComputeMEscnetModules(number = )` indexes `adj_matrices`:

| `number` | Matrix | Definition |
| --- | --- | --- |
| 1 | `pcor` | partial correlation only |
| 2 | `bayes_cor1` | `abs(pcor) * abs(cor)` |
| 3 | `bayes_cor2` | `pmax(pcor, 0) * pmax(cor, 0)` |
| 4 | `bayes_cor3` | min-max normalised `pcor * cor` |
| 5 | `bayes_cor4` | signed-hybrid normalised `pcor * cor` |
| 6 | `bayes_cor5` | z-scored `pcor * cor` |
| 7 | `bayes_cor6` | `((1 + pmax(pcor, 0)) / 2) * pmax(cor, 0)` |

## Benchmarking

```r
load("go.Rdata")   # data frame: column 2 = gene symbol, column 3 = GO term

out <- run_wgcna_analysis(seurat_obj, merged_data, cell_type = "INH", go = go)
out$summary_stats
out$plot
```

`run_wgcna_analysis()` builds the Mescnet network, an hdWGCNA network from the
same metacell matrix and a random network with the same genes, then scores all
three with EGAD neighbour-voting AUROC.

## Seurat v5 compatibility

`hdWGCNA::MetacellsByGroups()` and `hdWGCNA::MetaspotsByGroups()` predate Seurat
v5 and do not handle the new assay layers. `Mescnet` exports replacements with
the same names and signatures that

- read the expression matrix with `layer` on SeuratObject >= 5.0.0 and with
  `slot` otherwise,
- call `SeuratObject::JoinLayers()` after merging metacells or metaspots,
- write the `ident.group` identity into the metaspot metadata so that `Idents<-`
  succeeds, and
- keep the grouping vector and the metacell/metaspot list in sync when a group
  fails to aggregate.

Because the names collide with hdWGCNA, load `Mescnet` **after** `hdWGCNA`, or
call `Mescnet::MetacellsByGroups()` explicitly.

## Citation

If you use `Mescnet`, please cite the accompanying manuscript. See the
repository for the current reference.

## License

GPL (>= 3). The metacell and metaspot aggregation functions are derived from
hdWGCNA, which is released under the same license.
