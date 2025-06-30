library("RColorBrewer")
library("grDevices")

plot_heatmap_movida <- function(se,
                                features,
                                row_column = NULL,
                                pathway = NULL,
                                #  mydata, is it just the count matrix?
                                use_symbol = FALSE, # da togliere tutte le reference
                                cluster_rows = TRUE,
                                cluster_columns = FALSE,
                                center_mean = TRUE,
                                scale_row = FALSE,
                                winsorize_threshold = NULL,
                                plot_title = NULL,
                                export_data = FALSE,
                                ...) {
                                  
  # parameters check
  if (is.null(features) || length(features) == 0) {
    stop("features cannot be empty or NULL")
  }
  
  if (!is.null(winsorize_threshold)) {
    stopifnot(is.numeric(winsorize_threshold))
    stopifnot(winsorize_threshold >= 0)
  }

  # Select features to plot and check their presence
  if (!is.null(row_column)) {
    available_features <- features[features %in% rowData(se)[[row_column]]]
    feature_indices <- match(available_features, rowData(se)[[row_column]])
    heatmap_data <- assay(se)[feature_indices, , drop = FALSE]
  } else {
    available_features <- features[features %in% rownames(se)]
    heatmap_data <- assay(se)[available_features, , drop = FALSE]
  }

  # to avoid problems later, remove the ones non-expressed and with variance = 0
  to_remove <- apply(heatmap_data, 1, var) == 0
  heatmap_data <- heatmap_data[!to_remove, , drop = FALSE]

  if (nrow(heatmap_data) < 2) {
    warning("Creating a heatmap with only one gene...")
  }

  hm_name <- "Expression \nvalues"  

  # Handle pllotting of other data
  if (center_mean) {
    heatmap_data <- heatmap_data - rowMeans(heatmap_data)
    hm_name <- "Expression \nvalues"
  }

  if (scale_row) {
    heatmap_data <- t(scale(t(heatmap_data)))
    hm_name <- "Z-scores \nExpression \nvalues"
  }

  # If cutting extreme values
  if (!is.null(winsorize_threshold)) {
    # do the winsoring
    heatmap_data[heatmap_data < -winsorize_threshold] <- -winsorize_threshold
    heatmap_data[heatmap_data > winsorize_threshold] <- winsorize_threshold
  }

  # Generate title
  if (is.null(plot_title)) {
    title <- paste0("Signature heatmap - ", pathway)
  } else {
    title <- plot_title
  }

  # Export data if requested
  if (export_data) {
    return(heatmap_data)
  }

  if (nrow(heatmap_data) == 0 || ncol(heatmap_data) == 0) {
    warning("No data available for the heatmap. Please check your input features.")
    return(NULL)
  }

  # Generate heatmap
  ch <- ComplexHeatmap::Heatmap(
    matrix = heatmap_data,
    column_title = title,
    name = hm_name,
    col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
    rect_gp = grid::gpar(col = "white", lwd = 0.5),
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    row_labels = rownames(heatmap_data),
    # row_labels = ifelse(use_symbol && !is.null(rowData(se)$SYMBOL),
    #   rowData(se)$SYMBOL[match(rownames(heatmap_data), rownames(rowData(se)))],
    #   rownames(heatmap_data)
    # ),
    ...
  )

  return(ComplexHeatmap::draw(ch, merge_legend = TRUE))
  # return(ch)
}

# idea -> mettere solo la funzione per creare la heatamp, le altre operazioni nel overview.
