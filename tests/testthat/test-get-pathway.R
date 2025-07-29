test_that("get_pathway_features returns gene set from pathway for A_vs_B", {
  contrast <- "A_vs_B"
  if (contrast %in% model$get_contrasts()) {
    fea <- model$get_fea("transcriptomics", contrast)
    if (nrow(fea) > 0) {
      path <- rownames(fea)[1]
      genes <- model$get_pathway_features(pathway = path, contrast = contrast, source = "transcriptomics")
      expect_true(is.character(genes))
    } else {
      skip("No enrichment rows for A_vs_B")
    }
  } else {
    skip("Contrast A_vs_B not found")
  }
})

test_that("get_pathway_features returns NULL for non-existing pathway", {
  expect_null(model$get_pathway_features("non_existing_pathway", "A_vs_B", "proteomics"))
})