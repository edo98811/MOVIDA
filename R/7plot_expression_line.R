#' Plot Average Expression Line Plot for Selected Entities
#'
#' This function generates a line plot showing the average expression of selected entities (e.g., genes) across groups from a \code{SummarizedExperiment} object. It can also export the underlying plotting data if requested.
#'
#' @param entities A character vector of entity names (e.g., gene symbols) to plot.
#' @param se_object A \code{SummarizedExperiment} object containing the expression data and sample metadata.
#' @param group_var A character string specifying the column name in \code{colData(se_object)} to use for grouping samples. Default is \code{"group"}.
#' @param export_data Logical; if \code{TRUE}, the function returns the data used for plotting instead of the plot. Default is \code{FALSE}.
#' @param data_type A character string indicating the type of data (for annotation purposes). Default is \code{"unknown"}.
#'
#' @return A \code{ggplot} object showing the average expression per entity across groups, or a data frame with the plotting data if \code{export_data = TRUE}. If none of the entities are found, returns a plot with a warning message.
#'
#' @details
#' The function checks which entities are present in the count matrix of the provided \code{SummarizedExperiment} object. It calculates the average expression for each entity within each group and plots the results as lines. If any entities are missing, a warning is issued.
#'
#' @import ggplot2
#' @importFrom SummarizedExperiment assays colData
#' @export
plot_expression_line_movida <- function(entities, se_object, group_var = "group", export_data = FALSE, data_type = "unknown", mean_median = "mean") {

  if (is.null(entities) || length(entities) == 0) {
    return(ggplot() +
      annotate("text",
        x = 0.5, y = 0.5,
        label = "Select at least one feature to plot",
        size = 6, hjust = 0.5, vjust = 0.5
      ) +
      theme_void())
  }

  # Use count matrix from SummarizedExperiment
  count_matrix <- assays(se_object)$counts

  # Check which entities exist in count matrix
  existing_entities <- entities[entities %in% rownames(count_matrix)]
  missing_entities <- entities[!entities %in% rownames(count_matrix)]

  if (length(missing_entities) > 0) {
    warning(paste("Entities not found:", paste(missing_entities, collapse = ", ")))
  }

  if (length(existing_entities) == 0) {
    return(ggplot() +
      annotate("text",
        x = 0.5, y = 0.5,
        label = "No features found in data",
        size = 6, hjust = 0.5, vjust = 0.5
      ) +
      theme_void())
  }

  # Extract expression data for existing entities
  expr_data <- count_matrix[existing_entities, , drop = FALSE]

  # Get group information
  group_info <- colData(se_object)[[group_var]]

  # Calculate average expression per gene per group
  avg_expr <- data.frame(
    Entity = character(),
    Group = character(),
    Average_Expression = numeric(),
    stringsAsFactors = FALSE
  )

  for (entity in existing_entities) {
    for (group in unique(group_info)) {
      group_samples <- which(group_info == group)

      avg_val <- switch(mean_median,
        "mean" = mean(expr_data[entity, group_samples]),
        "median" = median(expr_data[entity, group_samples])
      )

      mean(expr_data[entity, group_samples])
      avg_expr <- rbind(avg_expr, data.frame(
        Entity = entity,
        Group = group,
        Average_Expression = avg_val
      ))
    }
  }

  # Export the plotting data if requested
  if (export_data) {
    return(avg_expr)
  }
  
  # Create line plot
  ggplot(avg_expr, aes(x = Group, y = Average_Expression, color = Entity, group = Entity)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 16),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.y = element_blank()
    ) +
    labs(
      title = "Average by Group",
      x = "Group",
      color = "Entity"
    )
}
