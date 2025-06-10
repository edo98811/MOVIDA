library(shiny)
source("functions/plot_expression.R")
source("modules/module_plot_container.R")

mod_sidebar_overview_ui <- function(id) {
  ns <- NS(id)
  uiOutput(ns("sidebar_ui_module")) # Placeholder for the sidebar UI
}
mod_sidebar_overview_server <- function(id, dashboard_elements, selected_row_source, movida_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Generate dynamic module IDs reactively and safely
    transcriptomics_id <- reactive({
      req(selected_row_source$selected, selected_row_source$source)
      paste0("destination_trans;source_", selected_row_source$source, ";selected_", selected_row_source$selected)
    })

    proteomics_id <- reactive({
      req(selected_row_source$selected, selected_row_source$source)
      paste0("destination_prot;source_", selected_row_source$source, ";selected_", selected_row_source$selected)
    })

    # Reactive plot functions
    plot_geneplot_transcriptomics <- reactive({
      req(selected_row_source$selected, selected_row_source$source)
      source <- selected_row_source$source
      selected <- selected_row_source$selected

      return(function(export_data = FALSE) {
        plot_expression_for_transcriptomics(
          selected,
          movida_data$get_expression("transcriptomics"),
          movida_data$get_anno_df("transcriptomics"),
          use_gene_name = (source == "transcriptomics"),
          export_data = export_data
        )
      })
    })

    plot_geneplot_proteomics <- reactive({
      req(selected_row_source$selected, selected_row_source$source)
      source <- selected_row_source$source
      selected <- selected_row_source$selected

      return(function(export_data = FALSE) {
        plot_expression_for_proteomics(
          selected,
          movida_data$get_expression("proteomics"),
          movida_data$get_anno_df("proteomics"),
          use_gene_name = (source == "transcriptomics"),
          export_data = export_data
        )
      })
    })

    plot_heatmap_transcriptomics <- reactive({
      req(selected_row_source$selected, selected_row_source$source)
      selected <- selected_row_source$selected

      return(function(export_data = FALSE) {
        plot_heatmap_for_transcriptomics(
          selected,
          movida_data$get_expression("transcriptomics"),
          movida_data$get_anno_df("transcriptomics"),
          export_data = export_data
        )
      })
    })

    plot_heatmap_proteomics <- reactive({
      req(selected_row_source$selected, selected_row_source$source)
      selected <- selected_row_source$selected

      return(function(export_data = FALSE) {
        plot_heatmap_for_proteomics(
          selected,
          movida_data$get_expression("proteomics"),
          movida_data$get_anno_df("proteomics"),
          export_data = export_data
        )
      })
    })

    # Dynamically render UI
    output$sidebar_ui_module <- renderUI({
      req(transcriptomics_id(), proteomics_id())

      tagList(
        h4("Related Omics Overview"),
        mod_plot_ui(
          id = ns("transcriptomics_plot"),
          title = "",
          show_export_data_btn = TRUE,
          show_export_plot_btn = TRUE
        ),
        mod_plot_ui(
          id = ns("proteomics_plot"),
          title = "",
          show_export_data_btn = TRUE,
          show_export_plot_btn = TRUE
        )
      )
    })

    # Connect server logic to the dynamically generated IDs
    observeEvent(transcriptomics_id(), {
      mod_plot_server("transcriptomics_plot", transcriptomics_id(), plot_geneplot_transcriptomics(), dashboard_elements) # Pass reactive plot function and id as not reactive
    })

    observeEvent(proteomics_id(), {
      mod_plot_server("proteomics_plot", proteomics_id(), plot_geneplot_proteomics(), dashboard_elements)
    })
  })
}
