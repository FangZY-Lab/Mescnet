test_that("the public API is exported", {
  funs <- c(
    "CreaterMEscnet",
    "CreaterMEscnet_data",
    "pick_power_scale_free_power",
    "ComputeMEscnetModules",
    "plot_metanet_auc",
    "run_wgcna_analysis",
    "MetacellsByGroups",
    "MetaspotsByGroups"
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
