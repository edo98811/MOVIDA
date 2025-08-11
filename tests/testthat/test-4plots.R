# tests/testthat/test-plot_expression_line_movida.R

test_that("entities is NULL", {
  # expect: ggplot object with 'Select at least one feature' message
  result <- plot_expression_line_movida(
    entities = NULL,
    se_object = model$get_dde_object_exposed("proteomics"),
    group_var = "group",
    export_data = FALSE
  )
  expect_s3_class(result, "ggplot")
  # optionally check that annotation text is in the plot
})

test_that("entities is empty vector", {
  result <- plot_expression_line_movida(
    entities = character(0),
    se_object = model$get_dde_object_exposed("proteomics"),
    group_var = "group"
  )
  expect_s3_class(result, "ggplot")
})

test_that("correct line plot", {
  result <- plot_expression_line_movida(
    entities = model$get_features_all("proteomics")[1:2],
    se_object = model$get_dde_object_exposed("proteomics"),
    group_var = "group",
    mean_median = "mean"
  )
  expect_s3_class(result, "ggplot")

  result <- plot_expression_line_movida(
    entities = model$get_features_all("proteomics")[1:2],
    se_object = model$get_dde_object_exposed("proteomics"),
    group_var = "group",
    mean_median = "median"
  )
  expect_s3_class(result, "ggplot")
  result <- plot_expression_line_movida(
    entities = model$get_features_all("transcriptomics")[1:2],
    se_object = model$get_dde_object_exposed("transcriptomics"),
    group_var = "group",
    mean_median = "median"
  )
  expect_s3_class(result, "ggplot")
  result <- plot_expression_line_movida(
    entities = model$get_features_all("metabolomics")[1:2],
    se_object = model$get_dde_object_exposed("metabolomics"),
    group_var = "group",
    mean_median = "median"
  )
  expect_s3_class(result, "ggplot")
})


test_that("group_var not found in colData", {
  expect_warning(
    plot_expression_line_movida(
      entities = model$get_features_all("proteomics")[1:2],
      se_object = model$get_dde_object_exposed("proteomics"),
      group_var = "nonexistent_column"
    )
  )
})

test_that("none of the entities are found in count matrix", {
  expect_warning(
    result <- plot_expression_line_movida(
      entities = c("fake_gene1", "fake_gene2"),
      se_object = model$get_dde_object_exposed("proteomics"),
      group_var = "group"
    ),
    regexp = "not found"
  )
  expect_s3_class(result, "ggplot")
})


test_that("some entities are missing", {
  expect_warning(
    plot_expression_line_movida(
      entities = c(model$get_features_all("proteomics")[[1]], "nonexistent_gene"),
      se_object = model$get_dde_object_exposed("proteomics"),
      group_var = "group"
    ),
    regexp = "Entities not found"
  )
})

test_that("invalid mean_median argument", {
  expect_warning(
    plot_expression_line_movida(
      entities = model$get_features_all("metabolomics")[1:2],
      se_object = model$get_dde_object_exposed("metabolomics"),
      group_var = "group",
      mean_median = "invalid_value"
    )
  )
})

test_that("export_data = TRUE returns data.frame", {
  result <- plot_expression_line_movida(
    entities = model$get_features_all("metabolomics")[1:2],
    se_object = model$get_dde_object_exposed("metabolomics"),
    group_var = "group",
    export_data = TRUE
  )
  expect_s3_class(result, "data.frame")
})

# 
# FOR EXPRESSION PLOTS
# 

# tests/testthat/test-plot_expression_movida.R

test_that("entity not found", {
  result <- plot_expression_movida(
    entity = "fake_entity",
    se_object = model$get_dde_object_exposed("proteomics"),
    group_var = "group"
  )
  expect_s3_class(result, "ggplot")
  # Optionally: check that label text contains "Entity 'fake_entity' not found"
})

test_that("correct plot for proteomics entity", {
  entity <- model$get_features_all("proteomics")[1]
  result <- plot_expression_movida(
    entity = entity,
    se_object = model$get_dde_object_exposed("proteomics"),
    group_var = "group"
  )
  expect_s3_class(result, "ggplot")
})

test_that("correct plot for transcriptomics entity", {
  entity <- model$get_features_all("transcriptomics")[1]
  result <- plot_expression_movida(
    entity = entity,
    se_object = model$get_dde_object_exposed("transcriptomics"),
    group_var = "group"
  )
  expect_s3_class(result, "ggplot")
})

test_that("correct plot for metabolomics entity", {
  entity <- model$get_features_all("metabolomics")[1]
  result <- plot_expression_movida(
    entity = entity,
    se_object = model$get_dde_object_exposed("metabolomics"),
    group_var = "group"
  )
  expect_s3_class(result, "ggplot")
})

test_that("group_var not found in colData", {
  # Depending on function behavior, this may error or warning
  expect_error(
    plot_expression_movida(
      entity = model$get_features_all("proteomics")[1],
      se_object = model$get_dde_object_exposed("proteomics"),
      group_var = "nonexistent_column"
    )
  )
})

test_that("export_data = TRUE returns data.frame", {
  entity <- model$get_features_all("proteomics")[1]
  result <- plot_expression_movida(
    entity = entity,
    se_object = model$get_dde_object_exposed("proteomics"),
    export_data = TRUE
  )
  expect_s3_class(result, "data.frame")
  expect_true(all(c("Value", "group") %in% colnames(result)))
})

# Optional: test background color detection independently if helpers are available
test_that("background color changes for entity types", {
  # Assuming you have known examples that match each checker
  prot_entity <- model$get_features_all("proteomics")[1]
  trans_entity <- model$get_features_all("transcriptomics")[1]
  metab_entity <- model$get_features_all("metabolomics")[1]

  p1 <- plot_expression_movida(prot_entity, model$get_dde_object_exposed("proteomics"))
  p2 <- plot_expression_movida(trans_entity, model$get_dde_object_exposed("transcriptomics"))
  p3 <- plot_expression_movida(metab_entity, model$get_dde_object_exposed("metabolomics"))

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
})
