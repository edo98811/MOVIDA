library(R6)
library(data.table)
library(gargoyle)


HeatmapContainer <- R6Class("HeatmapContainer",
  inherit = MainContentInterface,  # Inherit from MainContentInterface
  public = list(
    omics_type = NULL,
    results_type = NULL,
    FDR = NULL

    #' Initialize the class with parameters
    initialize = function(omics_type, results_type, FDR) {
      self$omics_type <- omics_type
      self$results_type <- results_type
      self$FDR <- FDR
    },

    #' Generate the main output data table
    main_output = function() {

      movida_data $get_de(type = self$omics_type)[
        FDR < self$FDR, 
        .(ID, log2FoldChange, pvalue, FDR)
      ]

      # Add a button column to trigger a gargoyle event
      input_data[, trigger_button := shiny::actionButton(
        inputId = paste0("trigger_", .I),
        label = "Trigger Event",
        onclick = paste0("gargoyle::trigger('event_", .I, "')")
      )]
      
      self$data_table <- input_data
      return(self$data_table)
    },

    #' Export the data table to a CSV file
    export_data = function() {
      if (is.null(self$data_table) || nrow(self$data_table) == 0) {
        stop("Data table is empty or not initialized.")
      }
      
      shiny::showModal(shiny::modalDialog(
        title = "Save Data Table",
        shiny::textInput("file_path", "Enter file path:", value = ""),
        shiny::actionButton("save_btn", "Save"),
        easyClose = TRUE
      ))
      
      shiny::observeEvent(shiny::input$save_btn, {
        file_path <- shiny::input$file_path
        if (file_path == "") {
          shiny::showNotification("File path cannot be empty.", type = "error")
        } else {
          fwrite(self$data_table, file = file_path)
          shiny::showNotification(sprintf("Data successfully exported to '%s'", file_path), type = "message")
          shiny::removeModal()
        }
      })
    }
  )
)