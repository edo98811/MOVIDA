library(shiny)
library(gargoyle)

source("modules/module_table_container.R")

# UI function for the module
mod_overview_ui <- function(id) {
  ns <- NS(id) # Namespace for the module

  navset_pill(
    nav_panel(
      title = "Differential Expression Overview",
      br(),
      div(
        h1("Transcriptomics"),
        mod_table_ui(ns("transcriptomics")),
        hr(),
        h1("Proteomics"),
        mod_table_ui(ns("proteomics"))
      )
    ),
    nav_panel(
      title = "Enrichment Results Overview",
      br(),
      div(
        h1("Transcriptomics"),
        mod_table_ui(ns("enrich_transcriptomics")),
        hr(),
        h1("Proteomics"),
        mod_table_ui(ns("enrich_proteomics"))
      )
    ),
    nav_spacer(),
    nav_item(div(
      style = "align-items: start; padding: 0 10px;",
      "Select Contrast: "
    )),
    nav_item(div(
      style = " align-items: end; align-items: center;  width: 100%;",
      selectInput(ns("contrast_selector"), NULL, choices = NULL)
    ))
  )
}

# Server function for the module
mod_overview_server <- function(id, dashboard_elements, row_to_select, movida_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    contrasts <- movida_data$get_contrasts()

    updateSelectInput(session, "contrast_selector",
      choices = contrasts,
      selected = contrasts[1]
    )

    contrast <- reactive({
      input$contrast_selector
    })

    # Update the contrast in the row_to_select reactive value
    observeEvent(contrast(), {
      row_to_select$contrast <- contrast()
      row_to_select$selected <- NULL
      row_to_select$source <- NULL
    })

    # Create reactive functions that close over the current contrast
    get_proteomics_table <- function() {
      movida_data$get_de("proteomics", contrast(), FDRpvalue = NULL, FDRadj = NULL)
    }

    # Create reactive functions that close over the current contrast (is it necessary?)
    get_transcriptomics_table <- function() {
      movida_data$get_de("transcriptomics", contrast(), FDRpvalue = NULL, FDRadj = NULL)
    }

    # Create reactive functions for enrichment tables
    get_enrich_proteomics_table <- function() {
      movida_data$get_enrichment("proteomics", contrast())
    }

    get_enrich_transcriptomics_table <- function() {
      movida_data$get_enrichment("transcriptomics", contrast())
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

    # Call mod_table_server for Proteomics
    mod_table_server(
      id = "enrich_transcriptomics",
      main_table_function = get_enrich_transcriptomics_table, # When the function is called the refresh is triggered
      selected_row = row_to_select
    )

    # Call mod_table_server for Transcriptomics
    mod_table_server(
      id = "enrich_proteomics",
      main_table_function = get_enrich_proteomics_table,
      selected_row = row_to_select
    )
  })
}
