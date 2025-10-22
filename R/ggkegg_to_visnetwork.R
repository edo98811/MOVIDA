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
ggkegg_to_igraph <- function(path_id, organism = "mmu", de_results = NULL) {
  # Validate pathway ID format
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
  nodes_df$KEGG <- vapply(nodes_df$name, remove_kegg_prefix_str, FUN.VALUE = character(1))

  # --- 2. Color and style nodes and edges ---
  # Color based on de resultsß
  nodes_df <- add_results_to_nodes(nodes_df, de_results)

  # Nodes
  nodes_df <- add_colors_to_nodes(nodes_df)
  nodes_df <- style_nodes(nodes_df)

  # Edges
  edges_df <- style_edges(edges_df)

  # --- 3. Scale node coordinates ---
  nodes_df <- scale_dimensions(nodes_df, factor = 2.5)

  # Create igraph object
  kegg_igraph <- igraph::graph_from_data_frame(edges_df, directed = FALSE, vertices = nodes_df)

  return(kegg_igraph)

  # --- 4. Define legend edges ---
  # legend_elements <- create_legend_dataframe(edges_df)

  # --- 5. Build  ---
  # return(
  #   visNetwork::visNetwork(nodes_df, edges_df, main = pathway_name) %>%
  #     visNetwork::visPhysics(enabled = FALSE) %>%
  #     visNetwork::visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, selectedBy = "group")
  #   # visNetwork::visLegend(addEdges = legend_elements$edges, useGroups = FALSE) %>% # , addNodes = legend_elements$nodes
  #   # visNetwork::visGroups(groupname = unique(nodes_df$group)) %>%
  #   # visNetwork::visLayout(randomSeed = 42)
  # )
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


style_nodes <- function(nodes_df) {
  nodes_df$shape <- "box"
  nodes_df$fixed <- TRUE # to avid moving nodes
  nodes_df$x <- x
  nodes_df$y <- y
  nodes_df$widthConstraint <- nodes$width
  nodes_df$heightConstraint <- nodes$height

  nodes_df$shape <- ifelse(nodes_df$type == "compound", "ellipse", "box")
  nodes_df$title <- paste0("<b>", nodes_df$label, "</b><br>KEGG ID(s): ", nodes_df$KEGG)

  return(nodes_df)
}


style_edges <- function(edges_df) {
  browser()
  edge_style_map <- list(
    compound = list(color = "black", dashes = FALSE, arrows = "to", label = ""),
    hidden_compound = list(color = "lightgray", dashes = FALSE, arrows = "to", label = ""),
    activation = list(color = "red", dashes = FALSE, arrows = "to", label = ""),
    inhibition = list(color = "blue", dashes = FALSE, arrows = "tee", label = ""),
    expression = list(color = "red", dashes = TRUE, arrows = "to", label = ""),
    repression = list(color = "blue", dashes = TRUE, arrows = "tee", label = ""),
    indirect_effect = list(color = "gray", dashes = TRUE, arrows = "to", label = ""),
    state_change = list(color = "gray", dashes = TRUE, arrows = "", label = ""),
    binding_association = list(color = "black", dashes = TRUE, arrows = "", label = ""),
    dissociation = list(color = "gray", dashes = TRUE, arrows = "to", label = ""),
    phosphorylation = list(color = "black", dashes = FALSE, arrows = "to", label = "+p"),
    dephosphorylation = list(color = "black", dashes = FALSE, arrows = "to", label = "-p"),
    glycosylation = list(color = "black", dashes = FALSE, arrows = "to", label = "+g"),
    ubiquitination = list(color = "black", dashes = FALSE, arrows = "to", label = "+u"),
    methylation = list(color = "black", dashes = FALSE, arrows = "to", label = "+m"),
    others_unknown = list(color = "black", dashes = TRUE, arrows = "to", label = "?")
  )

  edges_df$subtype[is.na(edges_df$subtype) || !(edges_subtype %in% names(edge_style_map))] <- "others_unknown"
  edges_df$subtype <- gsub("/", "_", edges_df$subtype)

  edges_df$color <- sapply(edges$subtype, function(x) edge_style_map[[x]]$color)
  edges_df$dashes <- sapply(edges$subtype, function(x) edge_style_map[[x]]$dashes)
  edges_df$arrows <- sapply(edges$subtype, function(x) edge_style_map[[x]]$arrows)
  edges_df$label <- sapply(edges$subtype, function(x) edge_style_map[[x]]$label)
}


#' Color nodes based on differential expression results.
#' @param nodes_df Data frame of nodes with a 'KEGG' column.
#' @param de_results_list Named list of differential expression results.
#'
#' @return nodes_df with an added 'color', 'source' and 'value' column.
add_results_to_nodes <- function(nodes_df, de_results_list) {
  # this function combines all results into a single data frame
  results_combined <- combine_results_in_dataframe(de_results_list)

  # Iterate through nodes_df and results_combined to add the results to the nodes
  nodes_df <- add_results_nodes(nodes_df, results_combined)

  return(nodes_df)
}

#' Add color palettes
#' @param nodes_df Data frame of nodes with 'value' and 'source' columns.
#' @return nodes_df with colored nodes based on their values.
#' @import RcolorBrewer
add_colors_to_nodes <- function(nodes_df) {

  palettes <- c(
    "RdGy",
    "RdBu",
    "PuOr",
    "PRGn",
    "PiYG",
    "BrBG"
  )

  # Get unique sources
  sources <- unique(na.omit(nodes_df$source))

  # Apply a color palette to each source
  for (source_index in seq_along(sources)) {
    palette <- palettes[[((source_index - 1) %% length(palettes)) + 1]]
    palette_colors <- RColorBrewer::brewer.pal(n = 11, name = palette)
    palette_ramp <- colorRampPalette(palette_colors)

    nodes_to_color <- nodes_df[nodes_df$source == sources[source_index], , drop = FALSE]

    if (nrow(nodes_to_color) > 1) {
      max_val <- max(nodes_to_color$value)
      min_val <- min(nodes_to_color$value)
    } else if (nrow(nodes_to_color) == 1) {
      max_val <- abs(nodes_to_color$value[1])
      min_val <- -abs(nodes_to_color$value[1])
    } else {
      next
    }

    # Compute color scaling around the center, maybe a better way to do this?
    range_val <- max(abs(c(max_val, min_val)))

    # Normalize values from -range_val to +range_val → [0,1]
    scaled_values <- (nodes_to_color$value + range_val) / (2 * range_val)
    scaled_values <- pmin(pmax(scaled_values, 0), 1) # clamp to [0, 1]

    # Assign colors based on scaled values
    nodes_to_color$color <- palette_ramp(100)[as.numeric(cut(scaled_values, breaks = 100, include.lowest = TRUE))]

    # Update main data frame
    nodes_df$color[nodes_df$source == sources[source_index]] <- nodes_to_color$color
  }

  return(nodes_df)
}

#' Download and cache KEGG KGML files.
#' @param pathway_id KEGG pathway ID (e.g., "hsa04110").
#' @param cache_dir Directory to store cached KGML files (default: "kgml_cache").
#' @return Path to the cached KGML file.
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
