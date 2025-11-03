#' Check the structure and validity of each entry in de_results
#'
#' @param de_entry An entry from the de_results list
#' @param name The name of the entry in the de_results list (for error messages)
#' 
#' @return Invisible TRUE if all checks pass, otherwise stops with an error
check_de_entry <- function(de_entry, name) {

  # --- Structure checks ---
  if (!is.list(de_entry) || !all(c("de_table", "value_column", "feature_column") %in% names(de_entry))) {
    warning(paste0("Each entry in de_results must be a list with elements: de_table, value_column, feature_column (problem in '", name, "')"))
    return(FALSE)
  }
  if (!inherits(de_entry$de_table, "data.frame")) {
    warning(paste0("de_results[['", name, "']]$de_table must be a data frame"))
    return(FALSE)
  }
  # Check that value_column is character and present in de_table
  if (!is.character(de_entry$value_column) || !(de_entry$value_column %in% colnames(de_entry$de_table))) {
    warning(paste0("de_results[['", name, "']]$value_column must be a column name in de_results[['", name, "']]$de_table"))
    return(FALSE)
  }
  # Check feature_column is character and present in de_table or is "rownames"
  if (!is.character(de_entry$feature_column) || !(de_entry$feature_column %in% c(colnames(de_entry$de_table), "rownames"))) {
    warning(paste0("de_results[['", name, "']]$feature_column must be a column name in de_results[['", name, "']]$de_table or 'rownames'"))
    return(FALSE)
  }

  # --- Row check ---

  # if feature_column is not in the colnames, if it is null, if it is not "rownames" then error
  if (!(de_entry$feature_column %in% colnames(de_entry$de_table)) &&
    de_entry$feature_column != "rownames") { warning(paste("Column", de_entry$feature_column, "not found in de_results: ", name))
    return(FALSE)
  }

  # --- Check that table is correct ---
  if (is.null(de_entry$de_table) || !(de_entry$value_column %in% colnames(de_entry$de_table))) {
    warning("Invalid de_table or value_column in de_results: ", name)
    return(FALSE)
  }

  return(TRUE)
}