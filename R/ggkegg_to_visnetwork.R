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
    if (!is.list(de_results) || is.null(names(de_results)) || any(names(de_results) == "")) {
      warning("de_results must be a named list or NULL.")
      return(NULL)
    }
    # Keep only valid entries
    valid_entries <- de_results[vapply(names(de_results), function(name) check_de_entry(de_results[[name]], name), logical(1))]
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
#' @return nodes_df with an added 'color', 'source' and 'value' column.
color_nodes <- function(nodes_df, de_results_list) {

  # Initializes the three columns in the nodes dataframe
  nodes_df[,c("value", "color", "source")] <- NA

  # get kegg ids from nodes, some nodes may have multiple kegg ids separated by ";"
  kegg_ids_graph <- expand_keggs(nodes_df)
 
  # For each entry in de_results_list, map values to KEGG ids in the graph. Returns a list of dataframes with columns (node_id, de_value) or NULL if no mapping
  kegg_values_list <- lapply(de_results_list, match_de_values)

  # Merge the single mappings dataframe for each entry in de_results_list, adds source info
  nodes_df <- add_values_to_nodes(nodes_df, kegg_values_list)

  # Add colors based on values
  nodes_df <- add_colors_to_nodes(nodes_df)

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
