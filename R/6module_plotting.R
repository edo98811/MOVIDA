
# DRAFT -------




# UI function for the module
mod_plotting_ui <- function(id) {
  ns <- NS(id) # Namespace for the module

  sources <- c(
    transcriptomics = "Transcriptomics",
    proteomics = "Proteomics",
    metabolomics = "Metabolomics"
  )

  # Source selector
  selectInput(
    ns("source_selector"),
    "Select Data Source",
    choices = sources,
    selected = "transcriptomics"
  )
  selectInput(
    ns("features_search_box"),
    "Select Features to Plot",
    choices = NULL,
    multiple = TRUE
  )
}

mod_plotting_server <- function(id, dashboard_elements, row_to_select, movida_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    sources <- c("transcriptomics", "proteomics", "metabolomics")
    
    # Update choices for the single features search box
    updateSelectInput(
      session,
      "features_search_box",
      choices = movida_data$get_all_features(input$source_selector)
    )


    # Save selected features



    plot_function <- function(export_data = FALSE) {
      se <- movida_data$get_values_all(source <- input$source_selector, return_se = TRUE)
      features <- input$features_search_box

      plot_expression_movida_line(
        se,
        features,
        export_data = export_data
      )
    }

    # Set up plot server logic for this feature
    mod_plot_server(
      id = paste0(source, "_gene_line_", gsub(":", "", pathway)),
      plot_id = plot_id(),
      main_plotting_function = plot_function,
      dashboard_elements = dashboard_elements
    )
  })
}
