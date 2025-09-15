#' @title Plot Heatmap for MOVIDA Signatures
#' @description
#' Generates a heatmap for selected features (e.g., genes) from a SummarizedExperiment object, with options for centering, scaling, winsorizing, clustering, and exporting the processed data.
#'
#' @param se A \code{SummarizedExperiment} object containing the expression data.
#' @param features A character vector of feature names (e.g., gene names or IDs) to include in the heatmap.
#' @param row_column Optional. Name of the column in \code{rowData(se)} to use for matching features (e.g., gene symbols). If \code{NULL}, rownames of \code{se} are used.
#' @param pathway Optional. Name of the pathway or signature for labeling the plot.
#' @param use_symbol Logical. Deprecated. If \code{TRUE}, attempts to use gene symbols for row labels (currently not used).
#' @param cluster_rows Logical. If \code{TRUE}, cluster the rows (features) in the heatmap.
#' @param cluster_columns Logical. If \code{TRUE}, cluster the columns (samples) in the heatmap.
#' @param center_mean Logical. If \code{TRUE}, center each row by subtracting its mean.
#' @param scale_row Logical. If \code{TRUE}, scale each row to have zero mean and unit variance (z-score).
#' @param winsorize_threshold Optional. Numeric value. If provided, expression values are limited to the range \code{[-winsorize_threshold, winsorize_threshold]}.
#' @param plot_title Optional. Custom title for the heatmap plot.
#' @param export_data Logical. If \code{TRUE}, returns the processed heatmap data matrix instead of plotting.
#' @param ... Additional arguments passed to the underlying \code{heatmap} function.
#'
#' @details
#' The function filters out features with zero variance and optionally centers, scales, and winsorizes the data. It supports both base R and ComplexHeatmap plotting (the latter is commented out). If no data is available after filtering, a warning is issued and \code{NULL} is returned.
#'
#' @return
#' If \code{export_data = TRUE}, returns a matrix of processed expression values. Otherwise, generates a heatmap plot and returns \code{NULL}.
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats heatmap
#'
#' @examples
#' \dontrun{
#' plot_heatmap_movida(se, features = c("GeneA", "GeneB", "GeneC"))
#' }
#'
#' @export
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
    availablefeatures <- features[features %in% rowData(se)[[row_column]]]
    feature_indices <- match(availablefeatures, rowData(se)[[row_column]])
    heatmap_data <- assay(se)[feature_indices, , drop = FALSE]
  } else {
    availablefeatures <- features[features %in% rownames(se)]
    heatmap_data <- assay(se)[availablefeatures, , drop = FALSE]
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

  # Generate heatmap using base R graphics
  heatmap(
    x = as.matrix(heatmap_data),
    main = title,
    col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
    Rowv = if(cluster_rows) NULL else NA,
    Colv = if(cluster_columns) NULL else NA,
    labRow = rownames(heatmap_data),
    margins = c(5, 10),
    cexRow = 0.8,
    cexCol = 0.8,
    ...
  )
}

# iDEA -> mettere solo la funzione per creare la heatamp, le altre operazioni nel overview.
