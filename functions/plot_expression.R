# Boxplot of logCPM values for each group
library(hrbrthemes)
library(viridis)
library(ggplot2)
library(dplyr)
library(ggforce)
library(ggrepel)
library(SummarizedExperiment)

plot_expression_for_proteomics <- function(entity, se_object, anno_df, use_gene_name = FALSE, count_matrix = NULL, export_data = FALSE) {
  # input protein must be in this case a uniprot id

  # Use correct plot name
  if (use_gene_name) {
    if (!entity %in% anno_df$ENSEMBL) {
      return(ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "Correspondence not found", size = 6, hjust = 0.5, vjust = 0.5) +
        theme_void())
    }

    protein <- anno_df[anno_df$ENSEMBL == entity, "UNIPROT"]
  } else {
    protein <- entity
  }

  if (is.null(count_matrix)) {
    count_matrix <- assays(se_object)$counts
  } else if (!is.matrix(count_matrix)) stop("count matrix must be either null or a matrix of counts")

  plotting_data <- t(count_matrix[protein, , drop = FALSE])

  rownames(plotting_data) <- colnames(count_matrix)
  colnames(plotting_data) <- "Value"

  plotting_data <- merge(plotting_data, colData(se_object)[c("group")], by = "row.names")
  plotting_data$group <- as.factor(plotting_data$group)

  # Export the plotting data if requested
  if (export_data) {
    return(plotting_data)
  }

  # Return plot
  return(plotting_data %>%
    ggplot(aes(x = group, y = Value, fill = group)) +
    # geom_boxplot(outliers = FALSE) +
    # geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
    ggforce::geom_sina(aes(color = group), size = 1.5) +
    ggrepel::geom_text_repel(aes(label = Row.names), size = 2) +
    geom_boxplot(alpha = 0.3) +
    # scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 21),
      panel.background = element_rect(fill = "#ffe9e9", color = NA)
    ) +
    ggtitle(protein) +
    xlab("") +
    ylab(protein))
}


plot_expression_for_transcriptomics <- function(entity, se_object, anno_df, use_gene_name = FALSE, count_matrix = NULL, export_data = FALSE) {
  # input protein must be in this case a ensembl id

  # Use correct plot name
  if (use_gene_name) {
    gene <- entity
  } else {
    if (!entity %in% anno_df$UNIPROT) {
      return(ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "Correspondence not found", size = 6, hjust = 0.5, vjust = 0.5) +
        theme_void())
    }

    gene <- anno_df[anno_df$UNIPROT == entity, "ENSEMBL"]
  }

  # Use correct plot name
  if (is.null(count_matrix)) {
    count_matrix <- assays(se_object)$counts
  } else if (!is.matrix(count_matrix)) stop("count matrix must be either null or a matrix of counts")

  plotting_data <- t(count_matrix[gene, , drop = FALSE])

  rownames(plotting_data) <- colnames(count_matrix)
  colnames(plotting_data) <- "Value"
  plotting_data <- merge(plotting_data, colData(se_object)[c("group")], by = "row.names")
  plotting_data$group <- as.factor(plotting_data$group)

  # Export the plotting data if requested
  if (export_data) {
    return(plotting_data)
  }

  # Return plot
  return(plotting_data %>%
    ggplot(aes(x = group, y = Value, fill = group)) +
    # geom_boxplot(outliers = FALSE) +
    # geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
    ggforce::geom_sina(aes(color = group), size = 1.5) +
    ggrepel::geom_text_repel(aes(label = Row.names), size = 2) +
    geom_boxplot(alpha = 0.3) +
    # scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 21),
      panel.background = element_rect(fill = "#e8f8fe", color = NA)
    ) +
    ggtitle(gene) +
    xlab("") +
    ylab(gene))
}


plot_expression <- function(entity, se_object, export_data = FALSE, data_type = "unknown", group_var = "group") {

  # Use count matrix from SummarizedExperiment if not provided

  count_matrix <- assays(se_object)$counts

  # Check if entity exists in count matrix
  if (!entity %in% rownames(count_matrix)) {
    return(ggplot() +
      annotate("text",
        x = 0.5, y = 0.5,
        label = paste0("Entity '", entity, "' not found in data"),
        size = 6, hjust = 0.5, vjust = 0.5
      ) +
      theme_void())
  }

  # Create plotting data
  plotting_data <- t(count_matrix[entity, , drop = FALSE])
  rownames(plotting_data) <- colnames(count_matrix)
  colnames(plotting_data) <- "Value"

  # Merge with group information
  plotting_data <- merge(plotting_data, colData(se_object)[c(group_var)], by = "row.names")
  plotting_data$group <- as.factor(plotting_data[[group_var]])

  # Set display name
  display_name <- entity

  # Determine background color based on data type
  if (data_type == "proteomics") {
    bg_color <- "#ffe9e9" # Light red for proteomics
  } else if (data_type == "transcriptomics") {
    bg_color <- "#e8f8fe" # Light blue for transcriptomics
  } else if (data_type == "metabolomics") {
    bg_color <- "#ffffdd" # Light yellow for metabolomics
  }

  # Export the plotting data if requested
  if (export_data) {
    return(plotting_data)
  }

  # Return plot
  return(plotting_data %>%
    ggplot(aes(x = group, y = Value, fill = group)) +
    ggforce::geom_sina(aes(color = group), size = 1.5) +
    ggrepel::geom_text_repel(aes(label = Row.names), size = 2) +
    geom_boxplot(alpha = 0.3) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 21),
      panel.background = element_rect(fill = bg_color, color = NA)
    ) +
    ggtitle(display_name) +
    xlab("") +
    ylab(display_name))
}
