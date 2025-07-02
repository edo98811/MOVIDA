plot_expression_movida_line <- function(entities, se_object, group_var = "group", export_data = FALSE, data_type = "unknown") {
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
        label = "No entities found in data",
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
      avg_val <- mean(expr_data[entity, group_samples])
      avg_expr <- rbind(avg_expr, data.frame(
        Entity = entity,
        Group = group,
        Average_Expression = avg_val
      ))
    }
  }

  # Export the plotting data if requested
  if (export_data) {
    return(plotting_data)
  }

  # Create line plot
  ggplot(avg_expr, aes(x = Group, y = Average_Expression, color = Entity, group = Entity)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 16),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      title = "Average by Group",
      x = "Group",
      y = "Average Expression",
      color = "Gene"
    )
}
