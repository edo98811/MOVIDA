library(shiny)
library(DT)

# UI function
mod_table_ui <- function(id, dashboard = FALSE, show_export_data_btn = FALSE, show_export_plot_btn = FALSE) {
  ns <- NS(id)

  btns <- list()
  # Add export data button if the flag is set
  if (show_export_data_btn) {
    btns <- append(btns, actionButton(ns("btn_export_data"), "Export Data"))
  }

  # Define the layout with buttons and plot output
  bslib::page_fillable(
    layout_columns(btns), # Render buttons dynamically
    DTOutput(ns("table")) # Render the plot/table output
  )
}

# Server function
mod_table_server <- function(id, main_table_function, selected_row, table_type = "default") {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Render the data table
    data <- reactive({
      result <- as.data.frame(main_table_function())
      if (!is.data.frame(result)) {
        warning("Data is not a data frame - returning NULL")
        return(NULL)
      }
      result
    })

    # Observe row selection in the data table
    observeEvent(input$table_rows_selected, {
      selected <- input$table_rows_selected
      if (length(selected)) {
        selected_row$selected <- rownames(data())[selected]
        selected_row$source <- id
      }
    })

    # Check if 'data' is a data.frame (or tibble)
    output$table <- DT::renderDataTable({
      if (is.null(data())) {
        DT::datatable(data.frame())
      } else if (table_type == "minimal") {
        DT::datatable(
          data(),
          options = list(
            deferRender = TRUE,
            dom = "t",
            scrollY = "400px",
            paging = FALSE,
            ordering = FALSE,
            headerCallback = JS("function(thead, data, start, end, display){ $(thead).remove(); }")
          ),
          rownames = FALSE,
          colnames = NULL,
          selection = "single", # Enable row selection,
          class = "compact stripe"
        )
      } else {
        DT::datatable(
          data(),
          selection = "single"
          # extensions = "Buttons",
          # options = list(
          #   dom = "Bfrtip",
          #   buttons = list(
          #     list(
          #       extend = "csv",
          #       text = "Export CSV",
          #       filename = paste0("data_", Sys.Date())
          #     )
          #   )
          # )
        )
      }
    })

    # # Observe the export data button click
    # observeEvent(input$btn_export_data, {
    #   showModal(modalDialog(
    #     title = "Export Data",
    #     tagList(
    #       # Input for file path
    #       textInput(session$ns("file_path"), "Enter file path:", value = paste0("data_", Sys.Date(), ".csv")),
    #       actionButton(session$ns("confirm_export"), "Confirm"), # Confirm export button
    #       modalButton("Cancel") # Cancel button
    #     )
    #   ))
    # })

    # # Handle the export confirmation
    # observeEvent(input$confirm_export, {
    #   req(input$file_path) # Ensure file path is provided
    #   removeModal() # Close the modal dialog
    #   tryCatch(
    #     {
    #       # Attempt to write the data to the specified file path
    #       write.csv(main_table_function(), file = input$file_path, row.names = FALSE)
    #       showNotification(paste("Data successfully exported to", input$file_path), type = "message")
    #     },
    #     error = function(e) {
    #       # Show error notification if export fails
    #       showNotification(paste("Error exporting data:", e$message), type = "error")
    #     }
    #   )
    # })

    # observeEvent(input$btn_export_plot, {
    #   showModal(modalDialog(
    #     title = "Export Plot",
    #     textInput(session$ns("plot_file_path"), "Enter file path (e.g., plot.png):", value = "plot.png"),
    #     footer = tagList(
    #       modalButton("Cancel"),
    #       actionButton(session$ns("confirm_plot_export"), "Export")
    #     )
    #   ))
    # })

    # observeEvent(input$confirm_plot_export, {
    #   req(input$plot_file_path)
    #   removeModal()

    #   tryCatch(
    #     {
    #       # Assuming main_plotting_function() returns a ggplot object
    #       ggsave(filename = input$plot_file_path, plot = main_table_function(), width = 8, height = 6)

    #       showNotification(paste("Plot successfully exported to", input$plot_file_path), type = "message")
    #     },
    #     error = function(e) {
    #       showNotification(paste("Error exporting plot:", e$message), type = "error")
    #     }
    #   )
    # })
  })
}
