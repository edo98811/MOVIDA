#' Transform a ggkegg graph to a visNetwork object.
#'
#' @param path_id KEGG pathway ID (e.g., "hsa:04110" or "04110").
#' @param organism KEGG organism code (e.g., "hsa" for human, "mmu" for mouse).
#' @param de_results Named list of differential expression results. Each entry should be a list with elements: de_table (data.frame), value_column (character), feature_column (character), threshold (numeric).
#' @return A visNetwork object representing the pathway with colored nodes based on differential expression results.
#'
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#'
#' @export
ggkegg_to_visnetwork <- function(path_id, organism = "mmu", de_results = NULL) {

  if (!is_valid_pathway(path_id)) {
    stop("Invalid KEGG pathway ID format. Must be like 'hsa:04110' or '04110'.")
  }

  # Validate each entry in de_results
  if (!is.null(de_results)) {
    # Check de_results is a named list with required structure
    if (!is.list(de_results) || length(de_results) == 0) {
      stop("de_results must be a non-empty list")
    }

    # Check names and structure of each entry
    invisible(
      lapply(names(de_results), function(name) {
        entry <- de_results[[name]]
        if (!is.list(entry) || !all(c("de_table", "value_column", "feature_column") %in% names(entry))) {
          stop(paste0("Each entry in de_results must be a list with elements: data, value_column, feature_column (problem in '", name, "')"))
        }
        if (!inherits(entry$de_table, "data.frame")) {
          stop(paste0("de_results[['", name, "']]$de_table must be a data frame"))
        }
        if (!is.character(entry$value_column) || !(entry$value_column %in% colnames(entry$de_table))) {
          stop(paste0("de_results[['", name, "']]$value_column must be a column name in de_results[['", name, "']]$de_table"))
        }
        if (!is.character(entry$feature_column) || !(entry$feature_column %in% c(colnames(entry$de_table), "rownames"))) {
          stop(paste0("de_results[['", name, "']]$feature_column must be a column name in de_results[['", name, "']]$de_table or 'rownames'"))
        }
      })
    )
  }


  if (!is.character(organism) || length(organism) != 1) {
    stop("organism must be a single character string")
  }

  organism <- to_organism_kegg(organism)

  # first call downloads, later calls use cached version
  kgml_file <- download_kgml(path_id)

  # If download failed, return NULL
  if (is.null(kgml_file)) {
    warning("Failed to download KGML file for pathway ID: ", path_id)
    return(NULL)
  }

  # Plot_title
  pathway_name <- paste0("(", path_id, ") ", get_pathway_name(path_id))

  # --- 1. Get nodes and edges data frames ---
  nodes_df <- parse_kgml_entries(kgml_file)
  edges_df <- parse_kgml_relations(kgml_file)

  # Get kegg ids from nodes names
  # format" type:number (e.g., "hsa:1234", "cpd:C00022", "path:map00010")
  nodes_df$KEGG <- vapply(nodes_df$name, ids_from_nodes_id, FUN.VALUE = character(1))

  # --- 2. Color and style nodes and edges ---
  # Color based on de resultsß
  nodes_df <- color_nodes(nodes_df, de_results)

  # Nodes
  nodes_df <- style_nodes(nodes_df, organism)

  # Edges
  edges_df <- style_edges(edges_df)

  # --- 3. Scale node coordinates ---
  nodes_df <- scale_dimensions(nodes_df, factor = 2.5)

  # --- 4. Define legend edges ---
  # legend_elements <- create_legend_dataframe(edges_df)

  # --- 5. Build network with legend ---
  return(
    visNetwork::visNetwork(nodes_df, edges_df, main = pathway_name) %>%
      visNetwork::visPhysics(enabled = FALSE) %>%
      visNetwork::visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, selectedBy = "group")
    # visNetwork::visLegend(addEdges = legend_elements$edges, useGroups = FALSE) %>% # , addNodes = legend_elements$nodes
    # visNetwork::visGroups(groupname = unique(nodes_df$group)) %>%
    # visNetwork::visLayout(randomSeed = 42)
  )
}


#' Scale node dimensions for better visualization.
#' @param nodes_df Data frame of nodes with x and y coordinates.
#' @param factor Scaling factor (default: 2).
#' 
#' @return nodes_df with scaled x and y coordinates.
scale_dimensions <- function(nodes_df, factor = 2) {
  # Scale x and y coordinates to make the graph look nicer
  nodes_df$x <- as.numeric(nodes_df$x) * factor
  nodes_df$y <- -as.numeric(nodes_df$y) * factor # Invert y-axis

  return(nodes_df)
}


style_nodes <- function(nodes_df, organism, delete_na = FALSE) {

  return(nodes_df)
}

style_edges <- function(edges_df) {

  return(edges_df)
}

#' Color nodes based on differential expression results.
#' @param nodes_df Data frame of nodes with a 'KEGG' column.
#' @param de_results_list Named list of differential expression results.
#' 
#' @return nodes_df with an added 'color' and 'value' column.
color_nodes <- function(nodes_df, de_results_list) {

  # If no de_results provided, return nodes_df unchanged
  if (is.null(de_results_list) || length(de_results_list) == 0) {
    message("No de_results provided, nodes will not be colored based on differential expression.")
    return(nodes_df)
  }

  # For each entry in de_results_list, color nodes accordingly
  colors_mappings <- lapply(names(de_results_list), function(name) {

    # --- Extract DE results ---
    entry <- de_results_list[[name]]
    de_results <- entry$de_table
    value_column <- entry$value_column
    feature_column <- entry$feature_column
    de_results <- de_results[!is.na(de_results[[feature_column]]), c(feature_column, value_column), drop = FALSE]

    # --- Row check ---
    if (feature_column == "rownames") {
      de_results$rownames <- rownames(de_results)
    } else if (!(feature_column %in% colnames(de_results))) {
      stop(paste("Column", feature_column, "not found in de_results"))
    }

    # --- Check that table is correct ---
    if (is.null(de_results) && !value_column %in% colnames(de_results)) {
      message("de_results not being mapped on graph, if provided check that column names are correct")
      return(NULL)
    }

    # --- get id list ---
    # It could be that some entries have multiple KEGG ids separated by ";", this handles that situation
    KEGG_in_graph <- sapply(nodes_df$KEGG, function(id) {
      match(get_single_ids(id), ids_from_kegg_mappings(de_results[[feature_column]]))
    }, USE.NAMES = TRUE)

    # --- Map keg ids to de value ---
    # Returns named vector of de values
    KEGG_in_de <- sapply(KEGG_in_graph, get_de_value, USE.NAMES = TRUE)
    if (all(is.na(KEGG_in_de))) {
      message("All values for value_column are NA. No features to map.")
      return(NULL)
    } else {
      return(KEGG_in_de)
    }
  })

  # --- Merge all mappings and assign colors ---
  de_mapping_merged <- merge_lists(colors_mapping)
  de_mapping_merged_no_repetitions <- handle_repetitions(merged)
  de_mapping_with_colors <- color_with_values(vals, palette = c("blue", "white", "red"), na_color = "white")

  # --- Add colors to nodes_df --- 
  merge(de_mapping_with_colors, nodes_df, by.x = "names", by.y = "KEGG", all.y = TRUE)

  return(nodes_df)
}

#' Download and cache KEGG KGML files.
#' @param pathway_id KEGG pathway ID (e.g., "hsa04110").
#' @param cache_dir Directory to store cached KGML files (default: "kgml_cache").
#' @return Path to the cached KGML file.
#' @importFrom httr GET content http_error
#' 
#' @importFrom xml2 read_xml
download_kgml <- function(pathway_id, cache_dir = "kgml_cache") {

  # Ensure cache dir exists
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }

  # Local file path
  file_path <- file.path(cache_dir, paste0(pathway_id, ".xml"))

  # If already cached, return path
  if (file.exists(file_path)) {
    message("Using cached file: ", file_path)
    return(file_path)
  }

  # KEGG REST API URL
  url <- paste0("https://rest.kegg.jp/get/", pathway_id, "/kgml")

  # Download KGML file
  kgml_content <- get_kgml(url)

  ifelse(!is.null(content), writeBin(kgml_content, file_path), return(NULL))
  message("Downloaded and cached: ", file_path)

  return(file_path)
}
