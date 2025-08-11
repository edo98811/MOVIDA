# DRAFT -------

# UI function for the module
mod_plotting_ui <- function(id) {
  ns <- NS(id) # Namespace for the module

  sources <- c(
    "transcriptomics",
    "proteomics",
    "metabolomics"
  )

  page_fillable(
    accordion(
      accordion_panel("Line Plot", uiOutput(ns("line_plot"))),
      accordion_panel("Heatmap", uiOutput(ns("heatma_plot"))),
      accordion_panel("PCA", uiOutput(ns("pca_plot"))),
      id = ns("plotting_accordion"), open = TRUE, fillable = TRUE, title = "Available plots"
    )
  )
}
mod_plotting_server <- function(id, dashboard_elements, bookmarked_elements, movida_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    sources <- c(
      "transcriptomics",
      "proteomics",
      "metabolomics"
    )
    #### UI FOR LINE PLOT ####
    output$line_plot <- renderUI({
      layout_column_wrap(
        width = 1 / 2,
        style = "margin-top: 20px; margin-bottom: 20px;",
        div(

          ## selector
          # - median / mean
          # - source
          # - features_to_plot
          # - group_by
          # - groups to show

          # Median/Mean selector buttons
          div(
            style = "margin-bottom: 10px;",
            actionButton(
              ns("median_btn"),
              "Median",
              class = "btn-default"
            ),
            actionButton(
              ns("mean_btn"),
              "Mean",
              class = "btn-primary"
            )
          ),

          # Source selector
          selectInput(
            ns("source_selector"),
            "Select Data Source",
            choices = sources,
            selected = "transcriptomics"
          ),
          # Group column selector
          uiOutput(ns("group_column_ui")),
          uiOutput(ns("subset_group_ui")),
          uiOutput(ns("features_search_box_ui"))
        ),
        div(
          uiOutput(ns("dynamic_line_plot"))
        )
      )
    })

    ##### SERVER LOGIC FOR LINE PLOT #####

    # Reactive value to store mean/median selection
    mean_median <- reactiveVal("median")

    observeEvent(input$median_btn, {
      mean_median("median")
      updateActionButton(
        session,
        "median_btn",
        class = "btn-primary"
      )
      updateActionButton(
        session,
        "mean_btn",
        class = "btn-default"
      )
    })

    observeEvent(input$mean_btn, {
      mean_median("mean")
      updateActionButton(
        session,
        "mean_btn",
        class = "btn-primary"
      )
      updateActionButton(
        session,
        "median_btn",
        class = "btn-default"
      )
    })

    output$group_column_ui <- renderUI({
      req(input$source_selector)
      selectInput(
        ns("group_column"),
        "Group by:",
        choices = colnames(movida_data$get_metadata(input$source_selector)),
        selected = "group"
      )
    })

    output$subset_group_ui <- renderUI({
      req(input$source_selector, input$group_column)
      selectInput(
        ns("subset_group"),
        "Select Subset",
        choices = unique(movida_data$get_metadata(input$source_selector)[[input$group_column]]),
        selected = NULL,
        multiple = TRUE
      )
    })

    output$features_search_box_ui <- renderUI({
      selectInput(
        ns("features_search_box"),
        "Select features to Plot",
        choices = movida_data$get_features_all(input$source_selector),
        multiple = TRUE,
        selected = NULL
      )
    })

    plot_function <- reactive({
      req(
        input$features_search_box,
        input$subset_group,
        input$source_selector,
        input$group_column,
        mean_median()
      )

      return(function(export_data = FALSE) {
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
      })
    })

    plot_id <- reactive({
      paste0(
        "source_", input$source_selector,
        ";selected_", input$group_column,
        ";selected_metadata_subset",  paste(input$subset_group, collapse = "_"),
        ";selected_metadata_column", input$group_column,
        ";mean_median_", mean_median(),
        ";features_", paste(input$features_search_box, collapse = "_")
      )
    })

    output$dynamic_line_plot <- renderUI({
      mod_plot_ui(
        id = ns("dynamic_line_plot"),
        title = "Line Plot",
        show_export_data_btn = TRUE,
        show_export_plot_btn = TRUE
      )
    })

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
