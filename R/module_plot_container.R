library(shiny)
source("functions/plot_expression.R")
source("utils/utils_export_functions.R")
library(shinyFiles)
# UI function
mod_plot_ui <- function(id, title = "Pinned Plot", show_export_data_btn = TRUE, show_export_plot_btn = TRUE, show_export_code_btn = TRUE) {
  ns <- NS(id)

  # btns <- list(
  #   uiOutput(ns("dashboard_buttons")),
  # )

  page_fillable(
    div(
      style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;",
      tags$h4(title, style = "margin: 0;"),
      uiOutput(ns("dashboard_buttons_small"))
    ),
    plotOutput(ns("plot")),
    div(
      class = "btn-group w-100 py-2", role = "group",
      # style = "margin-top: 10px; display: flex; gap: 5px;",
      if (show_export_data_btn) {
        downloadButton(ns("btn_export_data"), "Data", class = "btn-sm", style = "flex: 1; height: 35px;", title = "Export Data")
      },
      if (show_export_plot_btn) {
        page_fillable(
          downloadButton(ns("btn_export_plot_svg"), "SVG", class = "btn-sm", style = "flex: 1; height: 35px;", title = "Export SVG"),
          downloadButton(ns("btn_export_plot_png"), "PNG", class = "btn-sm", style = "flex: 1; height: 35px;", title = "Export PNG")
        )
      },
      if (show_export_code_btn) {
        downloadButton(ns("btn_export_code"), "Code", class = "btn-sm", style = "flex: 1; height: 35px;", title = "Export Code")
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
mod_plot_server <- function(id, plot_id, main_plotting_function, dashboard_elements = NULL) {
  moduleServer(id, function(input, output, session) {
    message("Module started for: ", plot_id)

    # Render the main plot using the provided plotting function
    output$plot <- renderPlot({
      main_plotting_function()
    })

    # Render UI for dashboard buttons based on whether the plot is in the dashboard
    # output$dashboard_buttons <- renderUI({
    #   if (!is.null(dashboard_elements) && !is.null(dashboard_elements$elements[[plot_id]])) {
    #     # If the plot is already in the dashboard, show a "Remove from Dashboard" button
    #     actionButton(session$ns("btn_remove_dashboard"), "Remove from Dashboard")
    #   } else {
    #     # Otherwise, show an "Add to Dashboard" button
    #     actionButton(session$ns("btn_dashboard"), "Add to Dashboard")
    #   }
    # })
    output$dashboard_buttons_small <- renderUI({
      if (!is.null(dashboard_elements) && !is.null(dashboard_elements$elements[[plot_id]])) {
        # If the plot is already in the dashboard, show a "Remove from Dashboard" button
        div(
          actionButton(session$ns("btn_remove_dashboard"), "-", class = "btn-outline-danger btn-sm", style = "width: 40px;"),
          style = "padding-bottom: 7px;"
        )
      } else {
        # Otherwise, show an "Add to Dashboard" button
        div(
          actionButton(session$ns("btn_add_dashboard"), "+", class = "btn-outline-success btn-sm", style = "width: 40px;"),
          style = "padding-bottom: 7px;"
        )
      }
    })

    # Add the plot to the dashboard when the "Add to Dashboard" button is clicked
    observeEvent(input$btn_add_dashboard, {
      req(main_plotting_function)
      if (!is.null(dashboard_elements)) {
        dashboard_elements$elements[[plot_id]] <- main_plotting_function
      }
    })

    # Remove the plot from the dashboard when the "Remove from Dashboard" button is clicked
    observeEvent(input$btn_remove_dashboard, {
      req(dashboard_elements)
      dashboard_elements$elements[[plot_id]] <- NULL
      dashboard_elements$elements <- dashboard_elements$elements[!sapply(dashboard_elements$elements, is.null)]
    })

    # Provide a download handler for exporting data as a CSV file
    output$btn_export_data <- downloadHandler(
      filename = function() {
        paste("exported_data-", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        req(main_plotting_function)
        # Call the plotting function with export_data = TRUE to get the data
        data <- main_plotting_function(export_data = TRUE)
        write.csv(data, file)
      }
    )

    # Provide a download handler for exporting the plot as a PNG file
    output$btn_export_plot_png <- downloadHandler(
      filename = function() {
        paste("exported_plot-", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
        req(main_plotting_function)
        # Save the plot to a PNG file
        png(file)
        print(main_plotting_function())
        dev.off()
      }
    )

    # Provide a download handler for exporting the plot as a PNG file
    output$btn_export_plot_svg <- downloadHandler(
      filename = function() {
        paste("exported_plot-", Sys.Date(), ".svg", sep = "")
      },
      content = function(file) {
        req(main_plotting_function)
        # Save the plot to a PNG file
        svg(file)
        print(main_plotting_function())
        dev.off()
      }
    )

    # Provide a download handler for exporting the code used to generate the plot
    output$btn_export_code <- downloadHandler(
      filename = function() {
        paste("exported_code-", Sys.Date(), ".R", sep = "")
      },
      content = function(file) {
        req(main_plotting_function)
        # Generate the R code for the plotting function
        func_text <- export_code(main_plotting_function, plot_id)
        writeLines(func_text, file)
      }
    )
  })
  # Create the hash of the function for the id.
  # # Convert the function to a string
  # func_string <- deparse(my_function)
  # # Collapse into a single string to avoid line breaks
  # func_text <- paste(func_string, collapse = "\n")
  # # Create a hash
  # hash <- digest(func_text, algo = "sha256")
}
