library("RColorBrewer")
library("grDevices")

gs_heatmap <- function(se,
                      #  mydata, is it just the count matrix?
                       res_enrich = NULL,
                       use_symbol = FALSE, # da togliere tutte le reference
                       geneset_id = NULL,
                       genelist = NULL,
                       FDR = 0.05,
                       cluster_rows = TRUE,
                       cluster_columns = FALSE,
                       center_mean = TRUE,
                       scale_row = FALSE,
                       winsorize_threshold = NULL,
                       anno_col_info = NULL,
                       plot_title = NULL,
                       ...) {
  if (!is.null(winsorize_threshold)) {
    stopifnot(is.numeric(winsorize_threshold))
    stopifnot(winsorize_threshold >= 0)
  }

  # check that the data would ideally be a DST, so that it is not the counts/normalized?

  # if(geneset in the results)
  #   pick the genes from there
  # else
  #   option to override the geneset by providing a list
  #
  # idea: multiselect with gene names - but in the UI
  # internal matching to the IDs (in this function we use the ids already)

  # rownames(res_enrich) <- res_enrich[["gs_id"]]
  if (!is.null(geneset_id)) {
    if (geneset_id %in% res_enrich[["gs_id"]]) {
      thisset_name <- res_enrich[geneset_id, "gs_description"]
      thisset_members <- unlist(strsplit(res_enrich[geneset_id, "gs_genes"], ","))
      thisset_members_ids <- rowData(se)$ENSEMBL[match(thisset_members, rowData(se)$SYMBOL)]
    }
  } else {
    # overridable via a list
    if (!all(genelist %in% rownames(se))) {
      not_there <- genelist[!(genelist %in% rownames(se))]
      warning(
        "Some of the provided gene ids were not found in the SummarizedExperiment",
        "\nNot found: ", not_there
      )
    }
    thisset_members_ids <- intersect(genelist, rownames(se))
    thisset_name <- "Custom list"
  }

  sig_to_keep <- (thisset_members_ids %in% rownames(se)) #
  thisset_members_ids_available <- thisset_members_ids[sig_to_keep]

  mydata_sig <- assay(se)[thisset_members_ids_available, , drop = FALSE]
  # mydata_sig <- mydata[thisset_members_ids_available, , drop = FALSE]

  # to avoid problems later, remove the ones non-expressed and with variance = 0
  to_remove <- apply(mydata_sig, 1, var) == 0
  mydata_sig <- mydata_sig[!to_remove, , drop = FALSE]

  if (nrow(mydata_sig) < 2) {
    warning("Creating a heatmap with only one gene...")
  }

  hm_name <- "Expression \nvalues"

  if (center_mean) {
    mydata_sig <- mydata_sig - rowMeans(mydata_sig)
    hm_name <- "Expression \nvalues"
  }

  if (scale_row) {
    mydata_sig <- t(scale(t(mydata_sig)))
    hm_name <- "Z-scores \nExpression \nvalues"
  }

  extreme_value <- max(abs(range(mydata_sig)))
  if (!is.null(winsorize_threshold)) {
    # do the winsoring
    mydata_sig[mydata_sig < -winsorize_threshold] <- -winsorize_threshold
    mydata_sig[mydata_sig > winsorize_threshold] <- winsorize_threshold
  }
  # dim(mydata_sig)

  if (is.null(plot_title)) {
    title <- paste0("Signature heatmap - ", thisset_name, " - ", geneset_id)
  } else {
    title <- plot_title
  }

  ### anno_col_info <- anno_col_info[anno_col_info %in% colnames(colData(se))]
  ### sample_decoration <- as.data.frame(colData(se))[, anno_col_info, drop = FALSE]

  # could there be a way to make this programmatically & clever?

  ## if only one column: SO
  # anno_col_vals <- colData(se)[,anno_col_info,drop = TRUE]
  #
  # ha_cols <- list(
  #   Annotation = structure(
  #     brewer.pal(length(unique(anno_col_vals)), "Set1"),
  #     names = unique(as.character(anno_col_vals))
  #   )
  # )
  # deco_ha <- HeatmapAnnotation(
  #   name = "eheh",
  #   Annotation = anno_col_vals,
  #   col = ha_cols
  # )
  ### deco_ha <- HeatmapAnnotation(df = sample_decoration)
  # if(returnData) {
  #
  # }

  # pheatmap(mydata_sig,
  #          # annotation_col = anno_colData,
  #          cluster_rows = cluster_rows, cluster_cols = cluster_cols,
  #          scale = ifelse(scale_row, "row", "none"), main = title,
  #          labels_row = annotation_obj[rownames(mydata_sig), ]$gene_name,
  #          annotation_col = sample_decoration)

  #
  if (is.null(anno_col_info)) {
    ch <- ComplexHeatmap::Heatmap(
      matrix = mydata_sig,
      column_title = title,
      name = hm_name,
      col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
      rect_gp = grid::gpar(col = "white", lwd = 0.5),
      cluster_rows = cluster_rows,
      cluster_columns = cluster_columns,
      row_labels = ifelse(use_symbol && !is.null(rowData(se)$SYMBOL),
        rowData(se)$SYMBOL[match(rownames(mydata_sig), rownames(rowData(se)))],
        rownames(rowData(se))
      ),
      ...
    )
  } else {
    anno_col_info <- anno_col_info[anno_col_info %in% colnames(colData(se))]
    sample_decoration <- as.data.frame(colData(se))[, anno_col_info, drop = FALSE]

    # handling the color with user defined values
    col_list <- vector(mode = "list", length = ncol(sample_decoration))

    for (i in seq_len(ncol(sample_decoration))) {
      these_decos <- sample_decoration[, i, drop = TRUE]

      if (is.numeric(these_decos)) {
        max_val <- max(these_decos, na.rm = TRUE)
        min_val <- min(these_decos, na.rm = TRUE)
        mid_val <- median(these_decos, na.rm = TRUE)

        if (max_val * min_val >= 0) {
          # sequential palette, all the same sign
          these_cols <- circlize::colorRamp2(
            c(min_val, max_val),
            c("grey95", "darkred")
          )
        } else {
          # need for a diverging palette
          these_cols <- circlize::colorRamp2(
            c(min_val, mid_val, max_val),
            c("blue", "grey95", "red")
          )
        }
      }

      if (is.character(these_decos) | is.factor(these_decos)) {
        if (is.character(these_decos)) {
          these_decos <- factor(these_decos)
        }
        these_cols <- colorspace::rainbow_hcl(n = length(levels(these_decos)))
        names(these_cols) <- levels(these_decos)
      }

      col_list[[i]] <- these_cols
    }

    names(col_list) <- colnames(sample_decoration)

    deco_ha <- ComplexHeatmap::HeatmapAnnotation(df = sample_decoration, col = col_list)

    ch <- ComplexHeatmap::Heatmap(
      matrix = mydata_sig,
      column_title = title,
      name = hm_name,
      col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
      rect_gp = grid::gpar(col = "white", lwd = 0.5),
      cluster_rows = cluster_rows,
      cluster_columns = cluster_columns,
      row_labels = ifelse(use_symbol && !is.null(rowData(se)$SYMBOL),
        rowData(se)$SYMBOL[match(rownames(mydata_sig), rownames(rowData(se)))],
        rownames(rowData(se))
      ),
      top_annotation = deco_ha,
      ...
    )
  }
  # is this gggplot?
  ComplexHeatmap::draw(ch, merge_legend = TRUE)
}

# idea -> mettere solo la funzione per creare la heatamp, le altre operazioni nel overview.
