test_that("the public API is exported", {
  funs <- c(
    "CreaterMEscnet",
    "CreaterMEscnet_data",
    "pick_power_scale_free_power",
    "ComputeMEscnetModules",
    "plot_metanet_auc",
    "run_wgcna_analysis",
    "MetacellsByGroups",
    "MetaspotsByGroups",
    "Plot_UMAP",
    "Plot_Module_Scores",
    "Plot_MEscnet_Network",
    "Plot_Hub_Heatmap",
    "Mescnet_colors"
  )
  for (f in funs) {
    expect_true(is.function(get(f, envir = asNamespace("Mescnet"))), info = f)
    expect_true(f %in% getNamespaceExports("Mescnet"), info = f)
  }
})

test_that("the magrittr pipe is re-exported", {
  expect_identical(1:3 %>% sum(), 6L)
})

test_that("MetacellsByGroups and MetaspotsByGroups expose the hdWGCNA signature", {
  expect_true("wgcna_name" %in% names(formals(MetacellsByGroups)))
  expect_true("min_spots" %in% names(formals(MetaspotsByGroups)))
})

test_that("Mescnet_colors validates its input and returns hexadecimal colours", {
  expect_identical(Mescnet_colors(0), character())
  expect_length(Mescnet_colors(3), 3)
  expect_length(Mescnet_colors(25), 25)
  expect_true(all(grepl("^#", Mescnet_colors(3))))
  expect_error(Mescnet_colors(-1), "non-negative")
  expect_error(Mescnet_colors(1.5), "non-negative")
})

test_that("plotting functions expose the documented hub-selection parameters", {
  expect_identical(formals(Plot_MEscnet_Network)$top_n, 5)
  expect_identical(formals(Plot_Hub_Heatmap)$top_n, 5)
  expect_true("target_cells" %in% names(formals(Plot_Hub_Heatmap)))
  expect_true("gene_list" %in% names(formals(CreaterMEscnet)))
})

test_that("hub genes are selected per module by absolute kME", {
  kme <- data.frame(
    Module1 = c(0.1, -0.9, 0.5, 0.2),
    Module2 = c(-0.8, 0.1, 0.3, 0.7),
    row.names = paste0("gene", 1:4)
  )
  modules <- list(c("gene1", "gene2", "gene3"), c("gene4"))

  hubs <- Mescnet:::.Mescnet_hub_genes(
    kme, modules, rownames(kme), top_n = 2
  )

  expect_identical(hubs, c("gene2", "gene3", "gene4"))
})

# Minimal S4 stand-in so the association-matrix code can be exercised without
# a full Seurat object.
setClass("MescnetTestObject", slots = c(misc = "list"), where = environment())

test_that("pick_power_scale_free_power builds the BayesCor matrices", {
  skip_if_not_installed("ppcor")

  set.seed(1)
  datExpr <- matrix(stats::rnorm(80 * 10), nrow = 80, ncol = 10)
  colnames(datExpr) <- paste0("gene", seq_len(10))

  obj <- new(
    "MescnetTestObject",
    misc = list(MEscnet = list(datExpr = datExpr))
  )

  obj <- pick_power_scale_free_power(
    obj,
    powerVector = 1:6,
    RsquaredCut = 0.5,
    plot_diagnostics = FALSE
  )

  payload <- obj@misc$MEscnet

  expect_equal(dim(payload$MEscnet_ppcor), c(10, 10))
  expect_true(all(paste0("bayes_cor", 1:6) %in% names(payload)))
  expect_equal(names(payload$result_tables)[1], "pcor")
  expect_length(payload$recommended_powers, 7)

  # BayesCor scores must be finite, symmetric and free of NA.
  for (nm in paste0("bayes_cor", 1:6)) {
    m <- payload[[nm]]
    expect_false(anyNA(m), info = nm)
    expect_equal(m, t(m), tolerance = 1e-8, info = nm)
  }

  # bayes_cor1 is abs(pcor) * abs(cor), so it is non-negative.
  expect_true(all(payload$bayes_cor1 >= 0))
  # bayes_cor2 zeroes out negative partial correlations.
  expect_true(all(payload$bayes_cor2 >= 0))
})

test_that("ComputeMEscnetModules validates the requested network index", {
  obj <- new("MescnetTestObject", misc = list(MEscnet = list()))
  expect_error(ComputeMEscnetModules(obj), "adj_matrices")
})
