ids_from_nodes_id <- function(kegg_ids) {
  # Remove prefix (e.g., "cpd:", "mmu:", "ko:", "path:")
  separated_elements <- strsplit(kegg_ids, " ")
  ids <- lapply(separated_elements, function(x) sub("^[a-z]+:", "", x))
  ids <- paste(unlist(ids), collapse = ";")

  return(ids)
}

#' Convert KEGG IDs with prefixes to IDs without prefixes
#' @param kegg_ids Character vector of KEGG IDs with prefixes (e.g., "cpd:C00022", "mmu:1234", "ko:K00001 ko:K00002")
#' @return List of character vectors with KEGG IDs without prefixes
ids_from_kegg_mappings <- function(kegg_ids) {
  # Remove prefix (e.g., "cpd:", "mmu:", "ko:", "path:")
  ids <- lapply(kegg_ids, function(x) sub("^[a-z]+:", "", x))
  return(ids)
}


get_single_ids <- function(ids_string) {
  ids <- unlist(strsplit(ids_string, ";"))
  return(ids)
}

#' Map normalized values to colors using a specified palette
#' @param vals_norm Numeric vector of normalized values (e.g., between -1 and 1)
#' @param palette Character vector of colors to use for the gradient (e.g., c("blue", "white", "red"))
#' @param na_color Color to use for NA
#' @return Named Character vector of colors (feature: color)
color_with_values <- function(vals_norm, palette, na_color = "white") {

}

get_de_value <- function(i) {
  if (all(is.na(i))) {
    return(NA)
  } else {
    return(de_results[[value_column]][i[!is.na(i)][1]])
  }
}

#' Merge a list of character vectors into a single data frame, handling repetitions.
#' @param list_of_char_vectors A list where each element is a named character vector that contains KEGG IDs as names and associated values from each given results dataframe.
#' @return A data frame with unique KEGG IDs and their corresponding values.
merge_values <- function(list_of_char_vectors) {
  # Extract names and values, removing prefixes
  ids_list <- unlist(lapply(list_of_char_vectors, names))
  values_list <- unlist(list_of_char_vectors, use.names = FALSE)

  merged <- data.frame(id = ids_list, value = values_list, stringsAsFactors = FALSE)

  return(merged)
}

handle_repetitions <- function(df, method = "first") {
  if (method == "first") {
    df_unique <- df[!duplicated(df$id), ]
  } else {
    stop("Only method = 'first' is currently supported.")
  }
  return(df_unique)
}

is_valid_pathway <- function(pathway_id) {
  # Check if pathway_id matches KEGG pathway formats: "hsa04110" or "04110"
  if (!is.character(pathway_id) || length(pathway_id) != 1) {
    return(FALSE)
  }
  grepl("^[a-z]{2,3}\\d{5}$", pathway_id) || grepl("^\\d{5}$", pathway_id)
}

to_organism_kegg <- function(organism) {
  if (missing(organism) || !nzchar(organism)) {
    stop("You must provide a valid KEGG organism code.")
  }

  # Handle common abbreviations
  organism_code <- switch(tolower(organism),
    "hs" = "hsa",
    "mm" = "mmu",
    "human" = "hsa",
    "mouse" = "mmu",
    "hsa" = "hsa",
    "mmu" = "mmu",
    stop(sprintf("Unknown organism abbreviation: %s", organism))
  )

  return(organism_code)
}

# Function to get KEGG compound name
get_kegg_name <- function(id) {
  entry <- KEGGREST::keggGet(id)[[1]]
  return(entry$NAME[1]) # first name
}
