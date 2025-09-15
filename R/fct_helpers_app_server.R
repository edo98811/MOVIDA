
safe_value <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(NA_character_)
  }
  x
}


is_bookmarked <- function(bookmarks, feature, source, contrast) {
  if (nrow(bookmarks) == 0) return(FALSE)
  
  nrow(bookmarks[
    ((is.na(bookmarks$feature) & is.na(feature)) | bookmarks$feature == feature) &
    ((is.na(bookmarks$source) & is.na(source)) | bookmarks$source == source) &
    ((is.na(bookmarks$contrast) & is.na(contrast)) | bookmarks$contrast == contrast),
  ]) > 0
}


add_bookmark <- function(bookmarks, feature, source, contrast) {
  new_bookmark <- data.frame(
    feature = feature,
    source = source,
    contrast = contrast,
    stringsAsFactors = FALSE
  )
  rbind(bookmarks, new_bookmark)
}

remove_bookmark <- function(bookmarks, feature, source, contrast) {
  bookmarks[!(
    ((is.na(bookmarks$feature) & is.na(feature)) | bookmarks$feature == feature) &
    ((is.na(bookmarks$source) & is.na(source)) | bookmarks$source == source) &
    ((is.na(bookmarks$contrast) & is.na(contrast)) | bookmarks$contrast == contrast)
  ), ]
}

bookmark_ui_info <- function(bookmarks, feature, source, contrast) {
  bookmarked <- is_bookmarked(bookmarks, feature, source, contrast)
  list(
    label = if (bookmarked) "Remove from Bookmarks" else "Add to Bookmarks",
    color = if (bookmarked) "red" else "green"
  )
}
