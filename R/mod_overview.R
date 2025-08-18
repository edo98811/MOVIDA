#' module_plotting_lineplot UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_overview_ui <- function(id) {
  ns <- NS(id) # Namespace for the module

  navset_pill(
    nav_panel(
      title = "Overview",
      uiOutput(ns("EntitiesOverview")), # Placeholder for the sidebar UI
      uiOutput(ns("is_filtering_button")) # Placeholder for the sidebar UI
    ),
    nav_panel(
      title = "Differential Expression",
      uiOutput(ns("DEResults"))
    ),
    nav_panel(
      title = "Enrichment Results",
      uiOutput(ns("EnrichmentResults"))
    )
  )
}
# UI Design:

# The renderUI blocks inside lapply are functional but could benefit from clearer separation of concerns. Consider extracting repetitive UI generation logic into helper functions.
# Performance:

# The lapply loops for generating UI and server logic are efficient but could benefit from caching or memoization if the data sources are large or computationally expensive.
# Server function for the module
    
#' module_plotting_lineplot Server Functions
#'
#' @noRd 
mod_overview_server <- function(id, dashboard_elements, row_to_select, movida_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
      ns <- session$ns

    contrast_selector_ui <- function(ns, choices, name) {
      div(
        style = "display: flex; align-items: center; padding: 10px 10px 0 10px;",
        div(style = "margin-right: 10px;", "Select Contrast:"),
        div(
          style = "flex-grow: 1;",
          selectInput(ns(paste0("contrast_selector_", name)), NULL, choices = choices, width = "100%")
        )
      )
    }
    
    # Reactive values to store
    contrast_de <- reactive({
      req(input$contrast_selector_de)
      input$contrast_selector_de
    })

    contrast_enrichment <- reactive({
      req(input$contrast_selector_enrichment)
      input$contrast_selector_enrichment
    })

    #  To have the two selectors syncornnized
    observeEvent(input$contrast_selector_de, ignoreInit = TRUE, {
      if (!identical(input$contrast_selector_enrichment, input$contrast_selector_de)) {
        updateSelectInput(session, "contrast_selector_enrichment",
          selected = input$contrast_selector_de
        )
      }
    })

    observeEvent(input$contrast_selector_enrichment, ignoreInit = TRUE, {
      if (!identical(input$contrast_selector_de, input$contrast_selector_enrichment)) {
        updateSelectInput(session, "contrast_selector_de",
          selected = input$contrast_selector_enrichment
        )
      }
    })

    is_filtering <- reactiveVal(FALSE)

    # Function to filter the tables if needed
    output$is_filtering_button <- renderUI({
      if (is.null(row_to_select$selected)) {
        actionButton(ns("is_filtering_button"), label = "Select row to filter")
      } else if (is_filtering()) {
        actionButton(ns("is_filtering_button"), label = paste("Showing only features related to (click to show all)", row_to_select$selected))
      } else {
        actionButton(ns("is_filtering_button"), label = paste("Click to show only features related to", row_to_select$selected))
      }
    })

    # Observe the is_filtering button click to toggle is_filtering
    observeEvent(input$is_filtering_button, {
      is_filtering(!is_filtering())
    })

    sources <- movida_data$get_sources()

    # Loop through sources to create server logic for DE and enrichment tables
    observeEvent(contrast_de(), {
      row_to_select$contrast <- contrast_de()
      row_to_select$selected <- NULL
      row_to_select$source <- NULL
    })

    get_tables <- function(source) {
      if (is_filtering()) {
        req(row_to_select$selected)
        data <- movida_data$get_related_features(row_to_select$selected, source)
      } else {
        data <- movida_data$get_features_all(source)
      }

      data <- data.frame(features = data, stringsAsFactors = FALSE)
      rownames(data) <- data[, 1]

      search_box <- input[[paste0("feature_", source, "_search_box")]]
      if (!is.null(search_box) && search_box != "") {
        data <- data[grep(search_box, rownames(data), ignore.case = TRUE), , drop = FALSE]
      }

      if (is.null(data) || nrow(data) == 0) {
        return(NULL)
      }

      return(data)
    }

    lapply(names(sources), function(source) {
      data_function <- reactive({
        tryCatch(
          {
            get_tables(source)
          },
          error = function(e) {
            warning(paste("Error in main_table_function:", e$message))
            return(NULL)
          }
        )
      })

      mod_table_server(
        id = paste0("de_", source),
        main_table_function = reactive({
          movida_data$getDEA(source, contrast_de(), FDRpvalue = NULL, FDRadj = NULL)
        }),
        selected_row = row_to_select
      )

      mod_table_server(
        id = paste0("enrich_", source),
        main_table_function = reactive({
          movida_data$getFEA(source, contrast_enrichment())
        }),
        selected_row = row_to_select
      )

      mod_table_server(
        id = paste0("feature_", source),
        main_table_function = data_function,
        selected_row = row_to_select,
        table_type = "minimal"
      )
    })

    output$EnrichmentResults <- renderUI({
      div(
        contrast_selector_ui(ns, movida_data$get_contrasts(), "enrichment"),
        br(),
        div(
          lapply(names(sources), function(source) {
            div(
              h1(sources[[source]]),
              mod_table_ui(ns(paste0("enrich_", source))),
              hr()
            )
          })
        )
      )
    })

    output$DEResults <- renderUI({
      div(
        contrast_selector_ui(ns, movida_data$get_contrasts(), "de"),
        br(),
        div(
          lapply(names(sources), function(source) {
            div(
              h1(sources[[source]]),
              mod_table_ui(ns(paste0("de_", source))),
              hr()
            )
          })
        )
      )
    })

    output$EntitiesOverview <- renderUI({
      do.call(layout_column_wrap, c(
        list(width = 1 / length(sources), heights_equal = "row"),
        lapply(names(sources), function(source) {
          div(
            br(),
            h4(sources[[source]]),
            div(
              style = "",
              textInput(ns(paste0("feature_", source, "_search_box")), label = "Search feature", placeholder = "Type to search...")
            ),
            mod_table_ui(ns(paste0("feature_", source)))
          )
        })
      ))
    })
  })
}
