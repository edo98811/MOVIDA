test_that("get_pathwayfeatures", {
  
  pathway <- rownames(model$getFEA("proteomics", "A_vs_B"))[[1]]
  result <- model$get_pathwayfeatures(pathway, "A_vs_B", "proteomics")

  expect_type(result, "character")
  expect_true(is.vector(result))
})

test_that("get_pathwayfeatures_not_existent", {

  pathway <- "nonexistent_pathway"
  suppressWarnings(result <- model$get_pathwayfeatures(pathway, "C_vs_B", "proteomics"))

  expect_true(is.null(result))
})

test_that("get_values_subset_metadata", {
  meta <- unique(model$get_metadata("proteomics")$group)[[1]]
  result <- model$get_values_subset_metadata(meta, "proteomics", column = 'group')
  result_se <- model$get_values_subset_metadata(meta, "proteomics", column = 'group', return_se = TRUE)

  expect_true(is.matrix(result))
  expect_s4_class(result_se, "SummarizedExperiment")
})

test_that("get_values_subset_metadata_not_existent", {
 
  expect_warning(result <- model$get_values_subset_metadata("C", "proteomics"))
  expect_true(is.null(result))
})

test_that("get_values_subset_metadata_empty", {

  # Test with null metadata
  meta <- NULL
  result <- model$get_values_subset_metadata(meta, "proteomics")
  expect_true(is.null(result))

  # Test with empty metadata as vector
  meta <- c()
  result_se <- model$get_values_subset_metadata(meta, "proteomics")
  expect_true(is.null(result_se))
})

test_that("get_values_samples", {
  meta <- model$get_metadata("proteomics")
  samples <- rownames(meta)[1:2]
  result <- model$get_values_samples(samples, "proteomics")
  dde <- model$get_values_samples(samples, "proteomics", return_se = TRUE)

  # Check if result is a matrix and has the correct rownames
  expect_equal(colnames(result), samples)
  expect_true(is.matrix(result))

  # Check when return_se is set to TRUE
  expect_s4_class(dde, "SummarizedExperiment")
  expect_equal(colnames(dde), samples)
  expect_true(all(colnames(dde) %in% rownames(meta)))
})

test_that("get_values_samples_not_existent", {
  samples <- c("nonexistent_sample1", "nonexistent_sample2")
  expect_warning(result <- model$get_values_samples(samples, "proteomics"))
  expect_warning(dde <- model$get_values_samples(samples, "proteomics", return_se = TRUE))

  expect_true(is.null(result))
  expect_true(is.null(dde)) 

  samples <- c(rownames(model$get_metadata("proteomics"))[1], "nonexistent_sample2")
  expect_warning(result <- model$get_values_samples(samples, "proteomics"))
  expect_warning(dde <- model$get_values_samples(samples, "proteomics", return_se = TRUE))

  # Only the real sample should be returned
  real_sample <- samples[1]
  expect_equal(colnames(result), real_sample)
  expect_true(is.matrix(result))

  expect_s4_class(dde, "SummarizedExperiment")
  expect_equal(colnames(dde), real_sample)

})

test_that("get_valuesfeatures", {
  features <- model$getfeatures_all("proteomics")[1:2]
  result <- model$get_values_features(features, "proteomics")
  dde <- model$get_values_features(features, "proteomics", return_se = TRUE)

  expect_equal(rownames(result), features)
  expect_true(is.matrix(result))

  expect_s4_class(dde, "SummarizedExperiment")
  expect_equal(rownames(dde), features)
})

test_that("get_valuesfeatures_not_existent", {
  features <- c("nonexistentfeature1", "nonexistentfeature2")
  expect_warning(result <- model$get_values_features(features, "proteomics"))
  expect_warning(dde <- model$get_values_features(features, "proteomics", return_se = TRUE))

  expect_true(is.null(dde))
  expect_true(is.null(result))

  features <- c(model$getfeatures_all("proteomics")[1], "nonexistentfeature2")
  expect_warning(result <- model$get_values_features(features, "proteomics"))
  expect_warning(dde <- model$get_values_features(features, "proteomics", return_se = TRUE))

  # Only the real feature should be returned
  realfeature <- features[1]
  expect_equal(rownames(result), realfeature)
  expect_true(is.matrix(result))

  expect_s4_class(dde, "SummarizedExperiment")
  expect_equal(rownames(dde), realfeature)
})

test_that("getDEA", {

  res <- model$getDEA("proteomics", "A_vs_B")
  expect_true(is.data.frame(res))
})

test_that("getDEA_not_existent", {

  suppressWarnings(res <- model$getDEA("proteomics", "ngaff"))
  expect_true(is.null(res))
})

test_that("getFEA", {

  res <- model$getFEA("proteomics", "A_vs_B")
  expect_true(is.data.frame(res))
})

test_that("getFEA_not_existent", {

  suppressWarnings(res <- model$getFEA("proteomics", "C_vs_B"))
  expect_true(is.null(res))
})
