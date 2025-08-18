library(testthat)
test_that("safe_value returns NA for NULL or empty", {
  expect_true(is.na(safe_value(NULL)))
  expect_true(is.na(safe_value(character(0))))
})

test_that("safe_value returns value if present", {
  expect_equal(safe_value("test"), "test")
  expect_equal(safe_value(1), 1)
})

test_that("is_bookmarked detects when bookmark exists", {
  bookmarks <- data.frame(
    feature = c("A", "B"),
    source = c("X", "Y"),
    contrast = c("C1", "C2"),
    stringsAsFactors = FALSE
  )
  
  expect_true(is_bookmarked(bookmarks, "A", "X", "C1"))
  expect_false(is_bookmarked(bookmarks, "A", "X", "C9"))
})

test_that("is_bookmarked works with NA matches", {
  bookmarks <- data.frame(
    feature = c(NA, "B"),
    source = c("X", NA),
    contrast = c("C1", "C2"),
    stringsAsFactors = FALSE
  )
  
  expect_true(is_bookmarked(bookmarks, NA, "X", "C1"))
  expect_true(is_bookmarked(bookmarks, "B", NA, "C2"))
})

test_that("add_bookmark adds new entries", {
  bookmarks <- data.frame(
    feature = "A",
    source = "X",
    contrast = "C1",
    stringsAsFactors = FALSE
  )
  
  updated <- add_bookmark(bookmarks, "B", "Y", "C2")
  expect_equal(nrow(updated), 2)
  expect_true(any(updated$feature == "B" & updated$source == "Y"))
})

test_that("remove_bookmark removes matching entries", {
  bookmarks <- data.frame(
    feature = c("A", "B"),
    source = c("X", "Y"),
    contrast = c("C1", "C2"),
    stringsAsFactors = FALSE
  )
  
  updated <- remove_bookmark(bookmarks, "A", "X", "C1")
  expect_equal(nrow(updated), 1)
  expect_false(any(updated$feature == "A"))
})

test_that("remove_bookmark works with NA values", {
  bookmarks <- data.frame(
    feature = c(NA, "B"),
    source = c("X", NA),
    contrast = c("C1", "C2"),
    stringsAsFactors = FALSE
  )
  
  updated <- remove_bookmark(bookmarks, NA, "X", "C1")
  expect_equal(nrow(updated), 1)
  expect_false(any(is.na(updated$feature) & updated$source == "X"))
})

test_that("bookmark_ui_info returns correct label and color", {
  bookmarks <- data.frame(
    feature = "A",
    source = "X",
    contrast = "C1",
    stringsAsFactors = FALSE
  )
  
  info1 <- bookmark_ui_info(bookmarks, "A", "X", "C1")
  expect_equal(info1$label, "Remove from Bookmarks")
  expect_equal(info1$color, "red")
  
  info2 <- bookmark_ui_info(bookmarks, "B", "Y", "C2")
  expect_equal(info2$label, "Add to Bookmarks")
  expect_equal(info2$color, "green")
})
