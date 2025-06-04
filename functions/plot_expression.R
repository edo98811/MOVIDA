# Boxplot of logCPM values for each group
library(hrbrthemes)
library(viridis)

plot_expression_for_proteomics <- function(entity, se_object, anno_df, use_gene_name = FALSE, count_matrix = NULL) {
  # input protein must be in this case a uniprot id

  # Use correct plot name
  if (use_gene_name) {

    if (!entity %in% anno_df$ENSEMBL_ID) {
      return(ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "Correspondence not found", size = 6, hjust = 0.5, vjust = 0.5) +
        theme_void())
    }

    protein <- anno_df[anno_df$ENSEMBL_ID == entity, "UNIPROT_ID"]
  } else {
    protein <- entity
  }

  if (is.null(count_matrix)) {
    count_matrix <- assays(se_object)$counts
  } else if (!is.matrix(count_matrix)) stop("count matrix must be either null or a matrix of counts")

  plotting_data <- t(count_matrix[protein,, drop = FALSE])

  rownames(plotting_data) <- colnames(count_matrix)
  colnames(plotting_data) <- "UNIPROT_ID"

  plotting_data <- merge(plotting_data, colData(se_object), by="row.names")
  plotting_data$group <- as.factor(plotting_data$group)
  
  return(plotting_data %>%
    ggplot(aes(x = group, y = UNIPROT_ID, fill = group)) +
    # geom_boxplot(outliers = FALSE) +
    # geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
    ggforce::geom_sina(aes(color = group), size = 1.5) + 
    ggrepel::geom_text_repel(aes(label = Row.names), size = 2) +
    geom_boxplot(alpha=0.3) +
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


plot_expression_for_transcriptomics <- function(entity, se_object, anno_df, use_gene_name = FALSE, count_matrix = NULL) {

  # input protein must be in this case a ensembl id

  # Use correct plot name
  if (use_gene_name) {
    gene <- entity
  } else {

    if (!entity %in% anno_df$UNIPROT_ID) {
      return(ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "Correspondence not found", size = 6, hjust = 0.5, vjust = 0.5) +
        theme_void())
    }

    gene <- anno_df[anno_df$UNIPROT_ID == entity, "ENSEMBL_ID"]
  }

  # Use correct plot name
  if (is.null(count_matrix)) {
    count_matrix <- assays(se_object)$counts
  } else if (!is.matrix(count_matrix)) stop("count matrix must be either null or a matrix of counts")

  plotting_data <- t(count_matrix[gene,, drop = FALSE])

  rownames(plotting_data) <- colnames(count_matrix)
  colnames(plotting_data) <- "ENSEMBL_ID"

  plotting_data <- merge(plotting_data, colData(se_object), by="row.names")
  plotting_data$group <- as.factor(plotting_data$group)
  
  return(plotting_data %>%
    ggplot(aes(x = group, y = ENSEMBL_ID, fill = group)) +
    # geom_boxplot(outliers = FALSE) +
    # geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
    ggforce::geom_sina(aes(color = group), size = 1.5) + 
    ggrepel::geom_text_repel(aes(label = Row.names), size = 2) +
    geom_boxplot(alpha=0.3) +
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