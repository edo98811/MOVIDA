library(shiny)
source("functions/plot_expression.R")
source("modules/module_plot_container.R")

mod_sidebar_overview_ui <- function(id) {
  ns <- NS(id)
  navset_pill(
    nav_panel(
      title = "Prot",
      div(
        div(class = "btn-group", role = "group",
          tags$button(type = "Select", class = "btn btn-primary", "Left"),
          tags$button(type = "Bookmark", class = "btn btn-primary", "Middle")
        ),
        uiOutput(ns("sidebar_ui_plots_prot")),
      )
    ),
    nav_panel(
      title = "Prot",
      div(
        div(class = "btn-group", role = "group",
          tags$button(type = "Select", class = "btn btn-primary", "Left"),
          tags$button(type = "Bookmark", class = "btn btn-primary", "Middle")
        ),
        uiOutput(ns("sidebar_ui_plots_prot")),
      )
    ),
    nav_panel(
      title = "Prot",
      div(
        div(class = "btn-group", role = "group",
          tags$button(type = "Select", class = "btn btn-primary", "Left"),
          tags$button(type = "Bookmark", class = "btn btn-primary", "Middle")
        ),
        uiOutput(ns("sidebar_ui_plots_prot")),
      )
    )
  )
}
mod_sidebar_overview_server <- function(id, dashboard_elements, selected_row_source, movida_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # a pill con 3 tabs
    # get related features always called for each navtab based on the selected feature
    # it ereturns either the same feature or all the related ones
    # i make a plot for each of the feature, with a button to add to bookmark and select the feature
    #
    #
    #
    #

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
