stest_that("get_contrasts includes A_vs_B", {
  result <- model$get_contrasts()
  expect_true("A_vs_B" %in% result)
})

test_that("get_dea returns a data.frame or NULL for A_vs_B", {
  contrasts <- model$get_contrasts()
  if ("A_vs_B" %in% contrasts) {
    res <- model$get_dea("proteomics", "A_vs_B")
    expect_true(is.data.frame(res) || is.null(res))
  } else {
    skip("Contrast A_vs_B not available")
  }
})

test_that("get_fea returns enrichment results for A_vs_B", {
  contrasts <- model$get_contrasts()
  if ("A_vs_B" %in% contrasts) {
    res <- model$get_fea("proteomics", "A_vs_B")
    expect_true(is.data.frame(res) || is.null(res))
  } else {
    skip("Contrast A_vs_B not available")
  }
})