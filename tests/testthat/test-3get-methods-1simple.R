test_that("get_contrasts", {
  contrasts <- model$get_contrasts()
  expect_true(is.character(contrasts))
  expect_equal(contrasts, c("A_vs_B", "B_vs_C"))
})

test_that("get_values_all", {
  m <- model$get_values_all("proteomics")
  expect_true(is.matrix(m))

  se <- model$get_values_all("proteomics", return_se = TRUE)
  expect_s4_class(se, "SummarizedExperiment")
})

test_that("get_features_all", {
  feat <- model$get_features_all("proteomics")
  expect_true(is.character(feat))
})

test_that("get_metadata_columns", {
  cols <- model$get_metadata_columns("proteomics")
  expect_true(is.character(cols))
})

test_that("get_metadata", {
  meta <- model$get_metadata("proteomics")
  expect_true(is.data.frame(meta))
  expect_true("group" %in% colnames(meta))
  expect_equal(sort(unique(meta$group)), c("A", "B"))
})

test_that("get_dde_object", {

  to_test <- c("proteomics", "transcriptomics", "metabolomics")
  for (source in to_test) {
    dde <- model$get_dde_object_exposed(source)
    expect_s4_class(dde, "DeeDeeExperiment")
    expect_true(all(rownames(dde) %in% model$get_features_all(source)))
    expect_true(all(colnames(dde) %in% model$get_metadata(source)$sample_id))
  } 
  expect_error(model$get_dde_object_exposed("nonexistent_source"), "get_dde_object: Invalid source type.")
})