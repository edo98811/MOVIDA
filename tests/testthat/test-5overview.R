library(shinytest2)
library(testthat)

test_that("Bookmark button works correctly", {

  # Launch the app (replace with your app path if needed)
  app <- AppDriver$new(app_dir = ".", name = "app_server_test", seed = 123)

  # Initially, no bookmarks
  expect_equal(app$get_value(output = "selectedfeature")$children[[2]]$children[[1]]$attribs$label, "Add to Bookmarks")

  # Simulate user clicking the bookmark button
  app$click("bookmark_toggle")

  # After click, the label should change
  expect_equal(app$get_value(output = "selectedfeature")$children[[2]]$children[[1]]$attribs$label, "Remove from Bookmarks")

  # Toggle again to remove
  app$click("bookmark_toggle")
  expect_equal(app$get_value(output = "selectedfeature")$children[[2]]$children[[1]]$attribs$label, "Add to Bookmarks")

  # Simulate changing selection
  app$set_inputs(selected_row_source_selected = "Feature2")
  app$set_inputs(selected_row_source_source = "SourceB")
  app$set_inputs(selected_row_source_contrast = "ContrastY")

  # Check that UI updates
  sidebar_ui <- app$get_value(output = "sidebar_ui")
  expect_true(grepl("Select a row to show information|Mock sidebar UI", as.character(sidebar_ui)))

  app$stop()
})

