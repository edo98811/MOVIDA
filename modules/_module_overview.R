library(shiny)
library(gargoyle)

source("modules/module_table_container.R")

# UI function for the module
mod_overview_ui <- function(id) {
  ns <- NS(id) # Namespace for the module
  div(
    selectInput(ns("contrast_selector"), "Select Contrast:", choices = NULL),

    h1("Transcriptomics"),
    mod_table_ui(ns("transcriptomics")),
    hr(),
    h1("Proteomics"),
    mod_table_ui(ns("proteomics"))
  )
}

# Server function for the module
mod_overview_server <- function(id, dashboard_elements, row_to_select, movida_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    contrasts <- movida_data$get_contrasts()

    updateSelectInput(session, "contrast_selector",
                      choices = contrasts,
                      selected = contrasts[1])

    contrast <- reactive({
      input$contrast_selector
    })

    # Update the contrast in the row_to_select reactive value
    observe({
      row_to_select$contrast <- contrast()
    })

    # Create reactive functions that close over the current contrast
    get_proteomics_table <- function() {
      movida_data$get_de("proteomics", contrast(), FDRpvalue = NULL, FDRadj = NULL)
    }

    # Create reactive functions that close over the current contrast (is it necessary?)
    get_transcriptomics_table <- function() {
      movida_data$get_de("transcriptomics", contrast(), FDRpvalue = NULL, FDRadj = NULL)
    }


    # Call mod_table_server for Proteomics
    mod_table_server(
      id = "transcriptomics",
      main_table_function = get_transcriptomics_table, # When the function is called the refresh is triggered
      selected_row = row_to_select
    )

    # Call mod_table_server for Transcriptomics
    mod_table_server(
      id = "proteomics",
      main_table_function = get_proteomics_table,
      selected_row = row_to_select
    )
  })
}
