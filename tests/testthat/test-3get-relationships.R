
test_that("get_inchi_to_ensembl", {
  result <- model$get_inchi_to_ensembl()
  expect_true(is.data.frame(result) || is.null(result))
})

test_that("get_inchi_to_uniprot", {
  result <- model$get_inchi_to_uniprot()
  expect_true(is.data.frame(result) || is.null(result))
})

test_that("get_uniprot_to_ensembl", {
  result <- model$get_uniprot_to_ensembl()
  expect_true(is.data.frame(result) || is.null(result))
})

test_that("get_related_features", {
  skip_if(is.null(model$get_inchi_to_uniprot()), "Skipping test for get_related_features as uniprot_to_inchi is not available.")
  
  feature <- rownames(model$get_features_all("proteomics"))[[1]]
  result <- model$get_related_features(feature, "metabolomics")
  
  expect_true(is.data.frame(result) || is.null(result))

  features <- rownames(model$get_features_all("proteomics"))[1:3]
  result <- model$get_related_features(features, "metabolomics")

  expect_true(is.data.frame(result) || is.null(result))
})

test_that("get_related_features_not_existent", {

  skip_if(is.null(model$get_inchi_to_uniprot()), "Skipping test for get_related_features as uniprot_to_inchi is not available.")
  feature <- "nonexistentfeature2"
  result <- model$get_related_features(feature, "metabolomics")
  
  expect_true(is.null(result))

  features <- c("nonexistentfeature1", "nonexistentfeature2")
  result <- model$get_related_features(features, "metabolomics")
  
  expect_true(is.null(result))
})
