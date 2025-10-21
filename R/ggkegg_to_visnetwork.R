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
  nodes_df <- color_nodes(nodes_df, de_results)

  # Nodes
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
  nodes_df$shape = "box"
  nodes_df$fixed = TRUE  # to avid moving nodes
  nodes_df$x = x
  nodes_df$y = y
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

  edges_df$subtype[is.na(edges_df$subtype)] <- "others_unknown"
  edges_df$subtype <- gsub("/", "_", edges_df$subtype )

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

  # get kegg ids from nodes, some nodes may have multiple kegg ids separated by ";"
  # kegg_ids_in_graph <- expand_results(nodes_df$KEGG) # tested 
  results_combined <- combine_results_dataframe(de_results_list) # tested

  nodes_df <- add_results_nodes(nodes_df, results_combined)

  # For each entry in de_results_list, map values to KEGG ids in the graph. Returns a list of dataframes with columns (KEGG, de_value or NULL if no mapping)
  kegg_ids_with_values <- lapply(de_results_list, match_de_values, kegg_ids_in_graph = kegg_ids_in_graph) # tested
  kegg_ids_with_values <- Reduce(function(x, y) merge(x, y, by = c("KEGG"), all = TRUE), kegg_ids_with_values)
  colnames(kegg_ids_with_values) <- c("KEGG", names(de_results_list))

  # # Warn if values are present in more than one column (excluding KEGG and value_to_plot)
  # value_counts <- rowSums(!is.na(kegg_ids_with_values[ , names(de_results_list)]))
  # if (any(value_counts > 1)) {
  #   warning("Some KEGG IDs have values in more than one differential expression source column.")
  # }

  # kegg_ids_with_values$value_to_plot <- rowMeans(kegg_ids_with_values[ , -1], na.rm = TRUE)


  # kegg_ids_with_values$source <- apply(kegg_ids_with_values[ , -c(1, ncol(kegg_ids_with_values))], 1, function(row) {
  #   sources <- names(de_results_list)[!is.na(row)]
  #   if (length(sources) == 0) {
  #     return(NA)
  #   } else {
  #     return(paste(sources, collapse = ";"))
  #   }
  # })
#  Example of kegg_ids_with_values:
#       KEGG    trans       prot      metabo value_to_plot source
# 1    C00076       NA         NA -0.04287046   -0.04287046 metabo
# 2    C00165       NA         NA          NA           NaN   <NA>
# 3    C00338       NA         NA          NA           NaN   <NA>
# 4    C00575 1.368602         NA          NA    1.36860228  trans
# 5    C01245       NA         NA          NA           NaN   <NA>
# 6  hsa04010       NA         NA          NA           NaN   <NA>
# 7  hsa04070       NA         NA          NA           NaN   <NA>
# an idea may be not to have value to plot and source, only create the sources column in the nodes_df. the rest is also added in the function so that i can delelte all the things that now are commented and put them in the add_values_to_nodes function' 
# idea for the function: iterate through nodes, see if each kegg id is present in the dataframe, if yes generate the sources string and add the value (maybe first if multiple?) column and sources column before
# to map to
#    id    name           type  link  reaction graphics_name label fgcolor bgcolor
#    <chr> <chr>          <chr> <chr> <chr>    <chr>         <chr> <chr>   <chr>  
#  1 19    cpd:C00338     comp… http… NA       C00338        C003… #000000 #FFFFFF
#  2 20    hsa:5923 hsa:… gene  http… NA       RASGRF1, CDC… RASG… #000000 #BFFFBF
#  3 21    hsa:11221 hsa… gene  http… NA       DUSP10, MKP-… DUSP… #000000 #BFFFBF
#  4 22    hsa:1845 hsa:… gene  http… NA       DUSP3, VHR... DUSP… #000000 #BFFFBF

  # Merge the single mappings dataframe for each entry in de_results_list, adds source info
  nodes_df <- add_values_to_nodes(nodes_df, kegg_ids_with_values)

  # Add colors based on values
  nodes_df <- add_colors_to_nodes(nodes_df)

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
