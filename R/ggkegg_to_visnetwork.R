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
ggkegg_to_igraph <- function(path_id, organism = "mmu", de_results = NULL, return_type = "igraph") {
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
    de_results <- de_results[vapply(names(de_results), function(name) check_de_entry(de_results[[name]], name), logical(1))]
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
  # Color based on de results
  results_combined <- combine_results_in_dataframe(de_results)
  nodes_df <- add_results_nodes(nodes_df, results_combined)

  # Nodes
  nodes_df <- add_colors_to_nodes(nodes_df)
  nodes_df <- style_nodes(nodes_df)

  # Edges
  edges_df <- style_edges(edges_df)

  # --- 3. Scale node coordinates ---
  nodes_df <- scale_dimensions(nodes_df, factor = 2.5)

  # --- 4. Create igraph object ---
  kegg_igraph <- igraph::graph_from_data_frame(edges_df, directed = FALSE, vertices = nodes_df)

  if (return_type == "igraph") {
    return(kegg_igraph)
  } else if (return_type == "visNetwork") {
    return(
      visNetwork::visIgraph(kegg_igraph, main = pathway_name) %>%
        visNetwork::visPhysics(enabled = FALSE)
    )
  }

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
  # Base visual settings
  nodes_df$shape <- ifelse(nodes_df$type == "compound", "ellipse", "box")
  nodes_df$fixed <- TRUE # prevent node movement
  nodes_df$node_name <- nodes_df$name
  nodes_df$name <- nodes_df$id

  nodes_df$widthConstraint <- as.numeric(nodes_df$width)
  nodes_df$heightConstraint <- as.numeric(nodes_df$height)

  nodes_df[!nodes_df$name == "undefined", ] # delete the groups (maybe later do something else with them)

  nodes_df$title <- paste0("<b>", nodes_df$label, "</b><br>KEGG ID(s): ", nodes_df$KEGG)

  return(nodes_df)
}


style_edges <- function(edges_df) {
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

  # https://builtin.com/data-science/and-in-r#:~:text=The%20single%20sign%20version%20%7C%20returns,first%20element%20of%20each%20vector.
  edges_df$subtype <- gsub("/", "_", edges_df$subtype)
  edges_df$subtype[is.na(edges_df$subtype) | !(edges_df$subtype %in% names(edge_style_map))] <- "others_unknown"

  edges_df$color <- vapply(edges_df$subtype, function(x) edge_style_map[[x]]$color, character(1))
  edges_df$dashes <- vapply(edges_df$subtype, function(x) edge_style_map[[x]]$dashes, logical(1))
  edges_df$arrows <- vapply(edges_df$subtype, function(x) edge_style_map[[x]]$arrows, character(1))
  edges_df$label <- vapply(edges_df$subtype, function(x) edge_style_map[[x]]$label, character(1))

  return(edges_df)
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


#' Add color palettes
#' @param nodes_df Data frame of nodes with 'value' and 'source' columns.
#' @return nodes_df with colored nodes based on their values.
#' @importFrom RColorBrewer brewer.pal
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
  valid_nodes <- nodes_df[!is.na(nodes_df$source), , drop = FALSE]

  # Apply a color palette to each source
  for (source_index in seq_along(sources)) {
    palette <- palettes[[((source_index - 1) %% length(palettes)) + 1]]
    palette_colors <- brewer.pal(n = 11, name = palette)
    palette_ramp <- colorRampPalette(palette_colors)
    nodes_to_color <- valid_nodes[valid_nodes$source == sources[source_index], , drop = FALSE]

    if (nrow(nodes_to_color) > 1) {
      range_val <- max(abs(nodes_to_color$value))
    } else if (nrow(nodes_to_color) == 1) {
      range_val <- abs(nodes_to_color$value[[1]])
    } else {
      next
    }
    # For ggplot: https://stackoverflow.com/questions/79132520/symmetric-colorbar-for-values-but-print-colorbar-for-actual-observed-values

    # Cut values using real numeric range (from -range_val to +range_val)
    # Generate breaks
    breaks_seq <- seq(-range_val, range_val, length.out = 101)

    # Use the cut function to assign colors based on the breaks (cut retrieves the index of the corresponding bin)
    # cut: https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/cut
    # may be useful to add general info in the dataframe: https://stackoverflow.com/questions/42217741/how-do-i-add-an-attribute-to-an-r-data-frame-while-im-making-it-with-a-function
    nodes_to_color$color <- palette_ramp(100)[
      as.numeric(cut(nodes_to_color$value, breaks = breaks_seq, include.lowest = TRUE))
    ]

    # Update main data frame
    valid_nodes$color[valid_nodes$source == sources[source_index]] <- nodes_to_color$color
  }
  nodes_df$color[!is.na(nodes_df$source)] <- valid_nodes$color

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
