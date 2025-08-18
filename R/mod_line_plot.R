
#' module_plotting_lineplot UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_line_plot_ui <- function(id) {
  ns <- NS(id)
  tagList(
    layout_column_wrap(
      width = 1/2,
      style = "margin-top: 20px; margin-bottom: 20px;",
      div(
        # Median/Mean selector buttons
        div(
          style = "margin-bottom: 10px;",
          actionButton(ns("median_btn"), "Median", class = "btn-default"),
          actionButton(ns("mean_btn"), "Mean", class = "btn-primary")
        ),

        # Source selector
        selectInput(ns("source_selector"), "Select Data Source", choices = NULL),

        # Dynamic UI placeholders
        uiOutput(ns("group_column_ui")),
        uiOutput(ns("subset_group_ui")),
        uiOutput(ns("features_search_box_ui"))
      ),
      div(
        mod_plot_ui(
          id = ns("dynamic_line_plot"),
          title = "Line Plot",
          show_export_data_btn = TRUE,
          show_export_plot_btn = TRUE
        )
      )
    )
  )
}
    
#' module_plotting_lineplot Server Functions
#'
#' @noRd 
mod_line_plot_server <- function(id, dashboard_elements, movida_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Get sources
    sources <- movida_data$get_sources()

    # Populate initial source choices
    observe({
      updateSelectInput(session, "source_selector",
                        choices = sources,
                        selected = "transcriptomics")
    })

    # Reactive value for mean/median toggle
    mean_median <- reactiveVal("median")

    observeEvent(input$median_btn, {
      mean_median("median")
      updateActionButton(session, "median_btn", class = "btn-primary")
      updateActionButton(session, "mean_btn", class = "btn-default")
    })

    observeEvent(input$mean_btn, {
      mean_median("mean")
      updateActionButton(session, "mean_btn", class = "btn-primary")
      updateActionButton(session, "median_btn", class = "btn-default")
    })

    # Group column selector
    output$group_column_ui <- renderUI({
      req(input$source_selector)
      selectInput(
        ns("group_column"),
        "Group by:",
        choices = colnames(movida_data$get_metadata(input$source_selector)),
        selected = "group"
      )
    })

    # Subset group selector
    output$subset_group_ui <- renderUI({
      req(input$source_selector, input$group_column)
      selectInput(
        ns("subset_group"),
        "Select Subset",
        choices = unique(movida_data$get_metadata(input$source_selector)[[input$group_column]]),
        multiple = TRUE
      )
    })

    # Features search box
    output$features_search_box_ui <- renderUI({
      req(input$source_selector)
      selectInput(
        ns("features_search_box"),
        "Select features to Plot",
        choices = movida_data$get_features_all(input$source_selector),
        multiple = TRUE
      )
    })

    # Plot function reactive
    plot_function <- reactive({
      req(
        input$features_search_box,
        input$subset_group,
        input$source_selector,
        input$group_column,
        mean_median()
      )

      function(export_data = FALSE) {
        features <- input$features_search_box
        se <- movida_data$get_values_subset_metadata(
          input$subset_group,
          input$source_selector,
          column = input$group_column,
          return_se = TRUE
        )

        plot_expression_line_movida(
          features,
          se,
          group_var = input$group_column,
          export_data = export_data,
          mean_median = mean_median()
        )
      }
    })

    # Plot ID (for bookmarking/dashboard)
    plot_id <- reactive({
      paste0(
        "source_", input$source_selector,
        ";selected_", input$group_column,
        ";selected_metadata_subset", paste(input$subset_group, collapse = "_"),
        ";selected_metadata_column", input$group_column,
        ";mean_median_", mean_median(),
        ";features_", paste(input$features_search_box, collapse = "_")
      )
    })

    # Launch plot server module
    observe({
      req(plot_id(), plot_function())
      mod_plot_server(
        id = "dynamic_line_plot",
        plot_id = plot_id(),
        main_plotting_function = plot_function(),
        dashboard_elements = dashboard_elements
      )
    })
  })
}
    
## To be copied in the UI
# mod_module_plotting_lineplot_ui("module_plotting_lineplot_1")
    
## To be copied in the server
# mod_module_plotting_lineplot_server("module_plotting_lineplot_1")
