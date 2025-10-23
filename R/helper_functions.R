remove_kegg_prefix_str <- function(kegg_ids) {
  # Remove prefix (e.g., "cpd:", "mmu:", "ko:", "path:")
  separated_elements <- strsplit(kegg_ids, " ")
  ids <- lapply(separated_elements, function(x) sub("^[a-z]+:", "", x))
  ids <- paste(unlist(ids), collapse = ";")
  return(ids)
}

#' Convert KEGG IDs with prefixes to IDs without prefixes
#' @param kegg_ids Character vector of KEGG IDs with prefixes (e.g., "cpd:C00022", "mmu:1234", "ko:K00001 ko:K00002")
#' @return List of character vectors with KEGG IDs without prefixes
remove_kegg_prefix <- function(kegg_ids) {
  # Remove prefix (e.g., "cpd:", "mmu:", "ko:", "path:")
  ids <- sapply(kegg_ids, function(x) sub("^[a-z]+:", "", x))
  return(ids)
}

#' Add results from combined results data frame to nodes data frame
#' @param nodes_df Data frame of nodes with a column 'KEGG'
#' @param results_combined Data frame with combined results containing columns: KEGG, value, source
#' @return Updated nodes data frame with added columns: value, color, source, text
add_results_nodes <- function(nodes_df, results_combined) {
  nodes_df[, c("value", "color", "source")] <- NA # Initialize new columns
  nodes_df$text <- ""

  # Iterate through nodes_df and results_combined to map values
  # For each node iterate over all results_combined
  for (i in seq_len(nrow(nodes_df))) {
    for (j in seq_len(nrow(results_combined))) {
      if (grepl(results_combined$KEGG[j], nodes_df$KEGG[i])) { # I made sure in both cases the keggs are without the prefix
        # If no value assigned to the node, assign the one from results_combined
        if (is.na(nodes_df$value[i])) {
          nodes_df$value[i] <- results_combined$value[j]
          nodes_df$source[i] <- results_combined$source[j]
          # nodes_df$text[i] <- list()
        } else { # If value warn
          warning(paste0("Multiple results mapped to node ", nodes_df$KEGG[i], ". Keeping the first occurrence to color the node."))
        }

        # This will be added in any case
        sep <- ","
        nodes_df$text[i] <- paste0(
          nodes_df$text[i],
          "Source: ", results_combined$source[j], sep,
          "Value: ", results_combined$value[j], sep,
          "Id: ", results_combined$KEGG[j], ";"
        )
      }
    }
  }

  return(nodes_df)
}

color_nodes <- function(nodes_df) {

}


#' Combine multiple differential expression results into a single data frame
#' @param results_list A named list where each element is a differential expression result containing a data frame (de_table), value column name (value_column), and feature column name (feature_column)
#' @return A combined data frame with columns: KEGG, value, source
combine_results_in_dataframe <- function(results_list) {
  results <- lapply(names(results_list), function(de_entry_name) {
    de_entry <- results_list[[de_entry_name]]
    de_table <- de_entry$de_table
    value_column <- de_entry$value_column
    feature_column <- de_entry$feature_column

    de_table[[feature_column]] <- remove_kegg_prefix(de_table[[feature_column]])

    data.frame(
      KEGG = de_table[[feature_column]],
      value = de_table[[value_column]],
      source = rep(de_entry_name, nrow(de_table)),
      stringsAsFactors = FALSE
    )
  })

  return(do.call(rbind, results))
}

#' Handle multiple KEGG IDs in a single string separated by ";"
#' @param kegg_string Character string of KEGG IDs separated by ";"
#' @return Data frame with a single column 'KEGG' containing individual KEGG IDs
expand_keggs <- function(kegg_df) {
  # Initialize empty vectors to store results
  ids_out <- c()
  kegg_out <- c()

  # Loop through each row of the data frame
  for (i in seq_len(nrow(kegg_df))) {
    # Split the KEGG string by ";"
    split_ids <- unlist(strsplit(kegg_df$KEGG[i], ";"))
    # Remove the prefix before ":" in each KEGG ID
    split_ids <- sub(".*:", "", split_ids)
    # Append the row ids and KEGG IDs
    ids_out <- c(ids_out, rep(kegg_df$id[i], length(split_ids)))
    kegg_out <- c(kegg_out, split_ids)
  }

  # Return the expanded data frame
  return(data.frame(id = ids_out, KEGG = kegg_out, stringsAsFactors = FALSE))
}

#' Merge a list of character vectors into a single data frame, handling repetitions.
#' @param list_of_char_vectors A list where each element is a named character vector that contains KEGG IDs as names and associated values from each given results dataframe.
#' @return A data frame with two columns: KEGG IDs and their corresponding values. (there can be repetitions)
add_values_to_nodes <- function(mapping_data_frame) {
  # Extract names and values, removing prefixes
  # Function body needed
  return(merged)
}

#' Handle repetitions in a data frame by keeping the first occurrence.
#' @param df Data frame with potential repetitions in the 'id' column.
#' @param method Method to handle repetitions. Currently only "first" is supported.
#' @return Data frame with repetitions removed based on the specified method.
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

#' Convert organism name or abbreviation to KEGG organism code
#' @param organism Organism name or abbreviation (e.g., "human", "hs", "mouse", "mm")
#' @return KEGG organism code (e.g., "hsa" for human, "mmu" for mouse)
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

#' Map differential expression values to KEGG IDs present in the graph
#' @param de_entry A single differential expression result entry containing de_table, value_column, and feature_column
#' @param kegg_ids_in_graph Character vector of KEGG IDs present in the graph
#' @return A data frame with columns: node_id (KEGG ID), value (mapped value)
match_de_values <- function(de_entry, kegg_ids_in_graph) {
  de_table <- de_entry$de_table
  value_column <- de_entry$value_column
  feature_column <- de_entry$feature_column

  de_table[[feature_column]] <- remove_kegg_prefix(de_table[[feature_column]])
  de_results <- data.frame(KEGG = de_table[[feature_column]], value = de_table[[value_column]], stringsAsFactors = FALSE)

  # The result of this function is a data frame with KEGG ids and a added column for values
  kegg_ids_in_graph_with_values <- merge(kegg_ids_in_graph, de_results, by = "KEGG", all.x = TRUE, all.y = FALSE)
  return(kegg_ids_in_graph_with_values)
}
