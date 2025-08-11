library(testthat)

# Assume the setup-tests.R already created objects like:
# movida_list (with dde_prot, dde_trans, dde_metabo), groups, contrasts, etc.

test_that("check_ensembl works correctly", {
  expect_true(check_ensembl(c("ENSG000001", "ENSMUSG000000")))
  expect_warning(expect_false(check_ensembl(c("ENSG000001", "ABC123"))))
})

test_that("check_goterm works correctly", {
  expect_true(check_goterm(c("GO:0008150", "GO:0003674")))
  expect_warning(expect_false(check_goterm(c("GO:0008150", "NOTAGO"))))
})

test_that("check_chebi works correctly", {
  expect_true(check_chebi(c("CHEBI:12345", "CHEBI:67890")))
  expect_warning(expect_false(check_chebi(c("CHEBI:12345", "12345"))))
})

test_that("check_uniprot works correctly", {
  expect_true(check_uniprot(c("P12345", "Q8N158")))
  expect_warning(expect_false(check_uniprot(c("P12345", "INVALID!"))))
})

test_that("check_inchi works correctly", {
  expect_true(check_inchi("ABCDEFGHIJKLMN-ABCDEFGHIJ-K"))
  expect_warning(expect_false(check_inchi("INVALIDINCHI")))
})

test_that("check_symbol", {
  expect_true(check_symbol(c("abc", "gene1")))
  expect_warning(expect_false(check_symbol(c("abc", "Gene-1"))))
})

test_that("check_organism", {
  expect_true(check_organism("Hs"))
  expect_true(check_organism("Mm"))
  expect_warning(expect_false(check_organism("Dm")))
})

# test_that("check_metadata", {
#   # Should not error if matching
#   expect_silent(check_metadata(movida_list_shared))

#   # Create mismatched metadata
#   ml_bad <- movida_list
#   rownames(ml_bad$movida_list_shared) <- paste0(rownames(ml_bad$movida_list_shared), "X")
#   expect_error(check_metadata(movida_list_shared))
# })

# test_that("check_movida_list", {
#   expect_true(check_movida_list(movida_list))

#   ml_bad <- movida_list
#   ml_bad$dde_prot <- NULL
#   ml_bad$dde_trans <- NULL
#   ml_bad$dde_metabo <- NULL
#   expect_error(check_movida_list(ml_bad))

#   ml_bad2 <- movida_list
#   ml_bad2$dde_prot <- data.frame(a=1)  # not SummarizedExperiment
#   expect_error(check_movida_list(ml_bad2))
# })

test_that("check_movida_list_dde", {
  # Assuming movida_list has valid dde objects
  expect_silent(check_movida_list_dde(movida_list))

  ml_bad <- movida_list
  ml_bad$dde_prot <- data.frame(a=1)
  expect_error(check_movida_list_dde(ml_bad))
})

test_that("check_contrast returns correct validation", {
  groups <- c("A", "B", "C")
  expect_true(check_contrast("A_vs_B", groups))
  expect_warning(result <- check_contrast("A_vs_D", groups))
  expect_false(result)
})

test_that("validate_contrasts returns only valid contrasts", {
  groups <- c("A", "B", "C")
  contrasts <- c("A_vs_B", "B_vs_C")
  expect_equal(validate_contrasts(contrasts, groups), c("A_vs_B", "B_vs_C"))
})
