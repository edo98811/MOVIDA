library(shiny)
source("functions/plot_expression.R")
library(shinyFiles)
# UI function
mod_plot_ui <- function(id, title = "Plot Module", show_export_data_btn = TRUE, show_export_plot_btn = TRUE, show_export_code_btn = TRUE) {
  ns <- NS(id)

  # btns <- list(
  #   uiOutput(ns("dashboard_buttons")),
  # )
  navset_card_tab(
    title = title,
    nav_panel(
      "Plot",
      plotOutput(ns("plot"))
    ),
    nav_panel(
      "Actions",
      uiOutput(ns("dashboard_buttons")),
      if (show_export_data_btn) {
        actionButton(ns("btn_export_data"), "Export Data")
      },
      if (show_export_plot_btn) {
        actionButton(ns("btn_export_plot"), "Export Plot")
      },
      if (show_export_code_btn) {
        actionButton(ns("btn_export_code"), "Export Code")
      }
    )
  )
  # if (show_export_data_btn) {
  #   btns <- append(btns, actionButton(ns("btn_export_data"), "Export Data"))
  # }

  # if (show_export_plot_btn) {
  #   btns <- append(btns, actionButton(ns("btn_export_plot"), "Export Plot"))
  # }

  # tagList(
  #     do.call(tagList, btns),

  # )
}

# Server function
mod_plot_server <- function(id, main_plotting_function, dashboard_elements = NULL) {
  moduleServer(id, function(input, output, session) {
    message("Module started: ", id)

    output$plot <- renderPlot({
      if (is.reactive(main_plotting_function)) {
        main_plotting_function()()
      } else {
        main_plotting_function()
      }
    })
    # Create the hash of the function for the id.
    # # Convert the function to a string
    # func_string <- deparse(my_function)
    # # Collapse into a single string to avoid line breaks
    # func_text <- paste(func_string, collapse = "\n")
    # # Create a hash
    # hash <- digest(func_text, algo = "sha256")
    output$dashboard_buttons <- renderUI({
      # req(main_plotting_function)
      if (!is.null(dashboard_elements) && !is.null(dashboard_elements$elements[[id]])) {
        actionButton(session$ns("btn_remove_dashboard"), "Remove from Dashboard")
      } else {
        actionButton(session$ns("btn_dashboard"), "Add to Dashboard")
      }
    })

    observeEvent(input$btn_dashboard, {
      req(main_plotting_function)
      if (!is.null(dashboard_elements)) {
        dashboard_elements$elements[[id]] <- main_plotting_function()
      }
    })

    observeEvent(input$btn_remove_dashboard, {
      req(dashboard_elements)
      dashboard_elements$elements[[id]] <- NULL
    })
  })
}
