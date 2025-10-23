test_that("expand_keggs handles multiple KEGG IDs correctly", {
  input <- data.frame(
    id = c(1, 2, 3),
    KEGG = c("hsa:1234;hsa:5678", "cpd:C00022", "ko:K00001;ko:K00002;ko:K00003"),
    stringsAsFactors = FALSE
  )
  expected_output <- data.frame(id = c(1, 1, 2, 3, 3, 3), KEGG = c("1234", "5678", "C00022", "K00001", "K00002", "K00003"))
  actual_output <- expand_keggs(input)
  expect_equal(actual_output, expected_output)
})


test_that("remove_kegg_prefix_str removes prefixes and handles multiple IDs", {
  input <- c("hsa:1234 hsa:5678", "cpd:C00022", "ko:K00001 ko:K00002 ko:K00003")
  expected_output <- c("1234;5678", "C00022", "K00001;K00002;K00003")
  actual_output <- vapply(input, remove_kegg_prefix_str, FUN.VALUE = character(1), USE.NAMES = FALSE)
  expect_equal(actual_output, expected_output)
})

test_that("parse_kgml_edges load relationsps correctly", {
  edges_df <- parse_kgml_relations(kgml_path)

  edges_df_expected <- tibble::as_tibble(read.csv(edges_df_path, sep = ";", colClasses = "character"))
  # sapply(colnames(edges_df), function(x) class(edges_df[[x]]))

  expect_equal(edges_df, edges_df_expected)
})

test_that("parse_kgml_entries load ndoes correctly", {
  nodes_df <- parse_kgml_entries(kgml_path)

  nodes_df_expected <- tibble::as_tibble(read.csv(nodes_df_path, sep = ";", colClasses = "character"))

  expect_equal(nodes_df, nodes_df_expected)
})

test_that("combine_results_in_dataframe correctly merges DE results", {
  # Run the function using your already-defined inputs
  result <- combine_results_in_dataframe(de_results_list)

  # Check that result is a data frame
  expect_true(is.data.frame(result))

  # Check expected columns
  expect_equal(colnames(result), c("KEGG", "value", "source"))

  # Check that 'source' column matches names of results_list
  expect_setequal(unique(result$source), names(de_results_list))

  # Check KEGG IDs no longer have the prefix (assuming remove_kegg_prefix removes 'path:')
  expect_false(any(grepl("^hsa:", result$KEGG)))

  # Check for non-missing values
  expect_false(any(is.na(result$KEGG)))
  expect_false(any(is.na(result$value)))

  # Check that all entries from each source are included
  expected_nrows <- sum(sapply(de_results_list, function(x) nrow(x$de_table)))
  expect_equal(nrow(result), expected_nrows)
})
