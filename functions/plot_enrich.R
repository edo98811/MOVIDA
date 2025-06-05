plot_enrichment_heatmap <- function(pathway, se_trans, annotation_df_trans, only_selected_contrast = TRUE, output) {

  # Filter the table based on the selected contrast (to change for deedeee)
  shaked_enrich_table <- GeneTonic::shake_topGOtableResult([[selected_contrast]]$topGO_tbl)

  if (selected_pathway %in% shaked_enrich_table$gs_id) {
    gene_list <- shaked_enrich_table[shaked_enrich_table$gs_id == selected_pathway, "gs_genes"]
    g_split <- strsplit(gene_list, split = ",", fixed = TRUE)[[1]]

    se <- reactive_values$se_trans

    # This will be eliminated when the input is standardized
    group_label_map <- c(
      "ATII" = "AT2",
      "Finerenon" = "finerenon",
      "SGLT2I" = "sglt2i",
      "ATII_Finerenon" = "AT2finerenon",
      "ATII_SGLT2I" = "AT2sglt2i",
      "ctrl" = "ctrl"
    )

    if (only_selected_contrast == TRUE) {
      s <- input$selected_contrast
      parts <- strsplit(s, "_vs_")[[1]]
      group1_label <- group_label_map[[parts[1]]]
      group2_label <- group_label_map[[parts[2]]]

      se <- se[, colData(se)$group %in% c(group1_label, group2_label)]
    }

    thisset_members_ids <- annotation_df_trans$ensembl_gene_id[match(g_split, annotation_df_trans$SYMBOL)]
    order <- rownames(colData(se)[order(colData(se)$group, decreasing = FALSE), ])
    order_numbers <- match(order, rownames(colData(se)))

    vst_data <- DESeq2::vst(assay(se))

    if (length(intersect(thisset_members_ids, rownames(vst_data)))) {
      message("creating plot")
    }

    return(gs_heatmap(
      mydata = vst_data,
      se = se,
      res_enrich = shaked_enrich_table,
      geneset_id = selected_pathway,
      annotation_obj = annotation_df_trans,
      scale_row = TRUE,
      show_row_dend = FALSE,
      show_column_names = TRUE,
      column_order = order_numbers,
      anno_col_info = "group"
    ))

  } else {
      return(ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "Correspondence not found", size = 6, hjust = 0.5, vjust = 0.5) +
        theme_void())
  }
}