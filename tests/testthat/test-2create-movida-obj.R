test_that("MovidaModel object is properly initialized", {
  skip("Skipping test for MovidaModel initialization as the model is not available.")
  movida_list_shared <- list(
    se_prot = se_prot,
    se_metabo = se_metabo,
    se_trans = se_trans,
    organism = "Mm",
    metadata = shared_metadata
  )
  expect_true(R6::is.R6(model))
})


test_that("MovidaModel object is properly initialized", {
  model <<- MovidaModel$new(movida_list)
  expect_true(R6::is.R6(model))
})

test_that("MovidaModel object is properly initialized", {
  model_se <<- expect_warning(MovidaModel$new(movida_list_se))
  expect_true(R6::is.R6(model_se))
})
