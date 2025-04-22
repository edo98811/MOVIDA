
library(igraph)
enrichment_map <- function(res_enrich,
                           gtl = NULL,
                           n_gs = 50,
                           gs_ids = NULL,
                           overlap_threshold = 0.1,
                           scale_edges_width = 200,
                           scale_nodes_size = 5,
                           color_by = "gs_pvalue") {
  if (!is.null(gtl)) {
    checkup_gtl(gtl)
    dds <- gtl$dds
    res_de <- gtl$res_de
    res_enrich <- gtl$res_enrich
    annotation_obj <- gtl$annotation_obj
  }

  if (!color_by %in% colnames(res_enrich)) {
    stop(
      "Your res_enrich object does not contain the ",
      color_by,
      " column.\n",
      "Compute this first or select another column to use for the color."
    )
  }

  n_gs <- min(n_gs, nrow(res_enrich))

  gs_to_use <- unique(
    c(
      res_enrich$gs_id[seq_len(n_gs)], # the ones from the top
      gs_ids[gs_ids %in% res_enrich$gs_id] # the ones specified from the custom list
    )
  )

  overlap_matrix <- GeneTonic::create_jaccard_matrix(res_enrich,
    n_gs = n_gs,
    gs_ids = gs_ids,
    return_sym = FALSE
  )

  rownames(overlap_matrix) <- colnames(overlap_matrix) <- res_enrich[rownames(overlap_matrix), "gs_description"]

  om_df <- as.data.frame(overlap_matrix)
  om_df$id <- rownames(om_df)

  omm <- tidyr::pivot_longer(om_df, seq_len(length(gs_to_use)))
  colnames(omm) <- c("gs_1", "gs_2", "value")
  # eliminate rows of diagonal...
  omm <- omm[omm$gs_1 != omm$gs_2, ]
  # ... and the ones from the other triangular portion
  omm <- omm[!is.na(omm$value), ]

  # omm <- reshape2::melt(overlap_matrix)
  # omm <- omm[omm$Var1 != omm$Var2, ]
  # omm <- omm[!is.na(omm$value), ]

  # use this to construct the graph
  emg <- graph_from_data_frame(omm[, c(1, 2)], directed = FALSE)

  E(emg)$width <- sqrt(omm$value * scale_edges_width)
  emg <- delete_edges(emg, E(emg)[omm$value < overlap_threshold])

  idx <- match(V(emg)$name, res_enrich$gs_description)

  gs_size <- res_enrich$gs_de_count[idx]

  V(emg)$size <- scale_nodes_size * sqrt(gs_size)
  V(emg)$original_size <- gs_size

  col_var <- res_enrich[idx, color_by]
  # the palette changes if it is z_score VS pvalue
  if (all(col_var <= 1) & all(col_var > 0)) { # likely p-values...
    col_var <- -log10(col_var)
    # V(g)$color <- colVar
    mypal <- (scales::alpha(
      colorRampPalette(RColorBrewer::brewer.pal(name = "YlOrRd", 9))(50), 0.8
    ))
    mypal_hover <- (scales::alpha(
      colorRampPalette(RColorBrewer::brewer.pal(name = "YlOrRd", 9))(50), 0.5
    ))
    mypal_select <- (scales::alpha(
      colorRampPalette(RColorBrewer::brewer.pal(name = "YlOrRd", 9))(50), 1
    ))
    
    V(emg)$color.background <- mosdef::map_to_color(col_var, mypal, symmetric = FALSE, 
                                         limits = range(na.omit(col_var)))
    V(emg)$color.highlight <- mosdef::map_to_color(col_var, mypal_select, symmetric = FALSE, 
                                        limits = range(na.omit(col_var)))
    V(emg)$color.hover <- mosdef::map_to_color(col_var, mypal_hover, symmetric = FALSE, 
                                    limits = range(na.omit(col_var)))
    
    V(emg)$color.background[is.na(V(emg)$color.background)] <- "lightgrey"
    V(emg)$color.highlight[is.na(V(emg)$color.highlight)] <- "lightgrey"
    V(emg)$color.hover[is.na(V(emg)$color.hover)] <- "lightgrey"
  } else {
    # e.g. using z_score or aggregated value
    if (prod(range(na.omit(col_var))) >= 0) {
      # gradient palette
      mypal <- (scales::alpha(
        colorRampPalette(RColorBrewer::brewer.pal(name = "Oranges", 9))(50), 0.8
      ))
      mypal_hover <- (scales::alpha(
        colorRampPalette(RColorBrewer::brewer.pal(name = "Oranges", 9))(50), 0.5
      ))
      mypal_select <- (scales::alpha(
        colorRampPalette(RColorBrewer::brewer.pal(name = "Oranges", 9))(50), 1
      ))
      
      V(emg)$color.background <- mosdef::map_to_color(col_var, mypal, symmetric = FALSE, 
                                           limits = range(na.omit(col_var)))
      V(emg)$color.highlight <- mosdef::map_to_color(col_var, mypal_select, symmetric = FALSE, 
                                          limits = range(na.omit(col_var)))
      V(emg)$color.hover <- mosdef::map_to_color(col_var, mypal_hover, symmetric = FALSE, 
                                      limits = range(na.omit(col_var)))
      V(emg)$color.background[is.na(V(emg)$color.background)] <- "lightgrey"
      V(emg)$color.highlight[is.na(V(emg)$color.highlight)] <- "lightgrey"
      V(emg)$color.hover[is.na(V(emg)$color.hover)] <- "lightgrey"
      
    } else {
      # divergent palette to be used
      mypal <- rev(scales::alpha(
        colorRampPalette(RColorBrewer::brewer.pal(name = "RdYlBu", 11))(50), 0.8
      ))
      mypal_hover <- rev(scales::alpha(
        colorRampPalette(RColorBrewer::brewer.pal(name = "RdYlBu", 11))(50), 0.5
      ))
      mypal_select <- rev(scales::alpha(
        colorRampPalette(RColorBrewer::brewer.pal(name = "RdYlBu", 11))(50), 1
      ))
      
      V(emg)$color.background <- mosdef::map_to_color(col_var, mypal, symmetric = TRUE, 
                                           limits = range(na.omit(col_var)))
      V(emg)$color.highlight <- mosdef::map_to_color(col_var, mypal_select, symmetric = TRUE, 
                                          limits = range(na.omit(col_var)))
      V(emg)$color.hover <- mosdef::map_to_color(col_var, mypal_hover, symmetric = TRUE, 
                                      limits = range(na.omit(col_var)))
      
      V(emg)$color.background[is.na(V(emg)$color.background)] <- "lightgrey"
      V(emg)$color.highlight[is.na(V(emg)$color.highlight)] <- "lightgrey"
      V(emg)$color.hover[is.na(V(emg)$color.hover)] <- "lightgrey"
    }
  }

  V(emg)$color.border <- "black"

  # additional specification of edge colors
  E(emg)$color <- "lightgrey"

  # re-sorting the vertices alphabetically
  rank_gs <- rank(V(emg)$name)
  emg <- permute(emg, rank_gs)

  return(emg)
}
