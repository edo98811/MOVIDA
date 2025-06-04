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
      paste0("plot_for_transcriptomics_from_", selected_row_source$source, "_", selected_row_source$selected)
    })

    proteomics_id <- reactive({
      req(selected_row_source$selected, selected_row_source$source)
      paste0("plot_for_proteomics_from_", selected_row_source$source, "_", selected_row_source$selected)
    })

    # Reactive plot functions
    plot_geneplot_transcriptomics <- reactive({
      req(selected_row_source$selected, selected_row_source$source)
      source <- selected_row_source$source
      selected <- selected_row_source$selected

      function() {
        plot_expression_for_transcriptomics(
          selected,
          movida_data$get_expression("transcriptomics"),
          movida_data$get_anno_df("transcriptomics"),
          use_gene_name = (source == "transcriptomics")
        )
      }
    })

    plot_geneplot_proteomics <- reactive({
      req(selected_row_source$selected, selected_row_source$source)
      source <- selected_row_source$source
      selected <- selected_row_source$selected

      function() {
        plot_expression_for_proteomics(
          selected,
          movida_data$get_expression("proteomics"),
          movida_data$get_anno_df("proteomics"),
          use_gene_name = (source == "transcriptomics")
        )
      }
    })

    # Dynamically render UI
    output$sidebar_ui_module <- renderUI({
      req(transcriptomics_id(), proteomics_id())

      tagList(
        h4("Related Omics Overview"),
        mod_plot_ui(
          id = ns(transcriptomics_id()),
          title = "Transcriptomics",
          show_export_data_btn = TRUE,
          show_export_plot_btn = TRUE
        ),
        mod_plot_ui(
          id = ns(proteomics_id()),
          title = "Proteomics",
          show_export_data_btn = TRUE,
          show_export_plot_btn = TRUE
        )
      )
    })

    # Connect server logic to the dynamically generated IDs
    observeEvent(transcriptomics_id(), {
      mod_plot_server(transcriptomics_id(), plot_geneplot_transcriptomics, dashboard_elements) # Pass reactive plot function and id as not rective
    })

    observeEvent(proteomics_id(), {
      mod_plot_server(proteomics_id(), plot_geneplot_proteomics, dashboard_elements)
    })
  })
}
