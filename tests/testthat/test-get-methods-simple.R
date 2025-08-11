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

test_that("getfeatures_all", {
  feat <- model$getfeatures_all("proteomics")
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
