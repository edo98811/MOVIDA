test_that("get_metadata_columns returns character vector", {
  cols <- model$get_metadata_columns("proteomics")
  expect_true(is.character(cols))
})

test_that("get_metadata returns data.frame", {
  meta <- model$get_metadata("proteomics")
  expect_true(is.data.frame(meta))
})
