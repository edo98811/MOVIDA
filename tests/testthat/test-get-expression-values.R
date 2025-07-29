test_that("get_values_all returns assay matrix or SE object", {
  m <- model$get_values_all("proteomics")
  expect_true(is.matrix(m))
  se <- model$get_values_all("proteomics", return_se = TRUE)
  expect_s4_class(se, "SummarizedExperiment")
})

test_that("get_values_features returns expected rows proteomics", {
  features <- rownames(model$get_values_all("proteomics"))[1:2]
  result <- model$get_values_features(features, "proteomics")
  expect_equal(rownames(result), features)
})

test_that("get_values_samples returns sample subset", {
  samples <- colnames(model$get_values_all("proteomics"))[1:2]
  result <- model$get_values_samples(samples, "proteomics")
  expect_equal(colnames(result), samples)
})

test_that("get_values_subset_metadata subsets correctly", {
  group_col <- "group"
  if (group_col %in% model$get_metadata_columns("proteomics")) {
    group_val <- as.character(colData(model$get_values_all("proteomics", TRUE))[[group_col]][1])
    se <- model$get_values_subset_metadata(subset = group_val, source = "proteomics")
    expect_s4_class(se, "SummarizedExperiment")
  } else {
    skip("No 'group' column found in metadata")
  }
})

test_that("get_features_list returns features", {
  feat <- model$get_features_list("proteomics")
  expect_true(is.character(feat))
})

test_that("get_features_list returns valid UniProt IDs for proteomics", {
  feats <- model$get_features_list("proteomics")
  expect_true(is.character(feats) || is.list(feats))
  expect_true(length(feats) > 0)
  expect_true(suppressWarnings(check_uniprot(feats)))
})

test_that("get_features_list returns valid Ensembl IDs or Symbol for transcriptomics", {
  feats <- model$get_features_list("transcriptomics")
  expect_true(is.character(feats) || is.list(feats))
  expect_true(length(feats) > 0)
  expect_true(suppressWarnings(check_ensembl(feats)) || check_symbol(feats))
})

test_that("get_features_list returns valid InChI IDs or Chebi IDs for metabolomics", {
  feats <- model$get_features_list("metabolomics")
  expect_true(is.character(feats) || is.list(feats))
  expect_true(length(feats) > 0)
  expect_true(suppressWarnings(check_inchi(feats) || check_chebi(feats)))
})
