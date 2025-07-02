
# UI function for the module
mod_overview_ui <- function(id) {
  ns <- NS(id) # Namespace for the module

  sources <- c(
    transcriptomics = "Transcriptomics",
    proteomics = "Proteomics",
    metabolomics = "Metabolomics"
  )

  contrast_selector_ui <- function(ns) {
    div(
      style = "display: flex; align-items: center; padding: 10px 10px 0 10px;",
      div(style = "margin-right: 10px;", "Select Contrast:"),
      div(style = "flex-grow: 1;", selectInput(ns("contrast_selector"), NULL, choices = NULL, width = "100%"))
    )
  }

  navset_pill(
    nav_panel(
      title = "Overview",
      # layout_column_wrap(
      #   width = 1 / 3,
      #   fill = FALSE,
      #   lapply(names(sources), function(source) {
      #     div(
      #       h4(sources[[source]]),
      #       textInput(ns(paste0("feature_", source, "_search_box")), label = "Search:", placeholder = "Type to search..."),
      #       mod_table_ui(ns(paste0("feature_", source)))
      #     )
      #   })
      # ),

      do.call(layout_column_wrap, c(
        list(width = 1 / 3, heights_equal = "row"),
        lapply(names(sources), function(source) {
          div(
            br(),
            h4(sources[[source]]),
            div(style = "",
              textInput(ns(paste0("feature_", source, "_search_box")), label = "Search feature", placeholder = "Type to search...")
            ),
            mod_table_ui(ns(paste0("feature_", source)))
          )
        })
      )),
      uiOutput(ns("is_filtering_button")) # Placeholder for the sidebar UI
    ),
    nav_panel(
      title = "Differential Expression",
      div(
        contrast_selector_ui(ns),
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
    ),
    nav_panel(
      title = "Enrichment Results",
      div(
        contrast_selector_ui(ns),
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
    )
  )
}
# UI Design:

# The renderUI blocks inside lapply are functional but could benefit from clearer separation of concerns. Consider extracting repetitive UI generation logic into helper functions.
# Performance:

# The lapply loops for generating UI and server logic are efficient but could benefit from caching or memoization if the data sources are large or computationally expensive.
# Server function for the module
mod_overview_server <- function(id, dashboard_elements, row_to_select, movida_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    contrasts <- movida_data$get_contrasts()

    if (length(contrasts) > 0) {
      updateSelectInput(session, "contrast_selector",
        choices = contrasts,
        selected = contrasts[1]
      )
    }

    contrast <- reactive({
      req(input$contrast_selector)
      input$contrast_selector
    })

    is_filtering <- reactiveVal(FALSE)

    # Update the contrast in the row_to_select reactive value
    observeEvent(contrast(), {
      row_to_select$contrast <- contrast()
      row_to_select$selected <- NULL
      row_to_select$source <- NULL
    })

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

    sources <- c("transcriptomics", "proteomics", "metabolomics")

    # Freeze the row_to_select$selected value when is_filtering is TRUE
    # If is_filtering is TRUE, the selected row will not change and the tables are not reactive aymore
    # observeEvent(is_filtering(), {
    #   if (is_filtering()) {
    #     freezeReactiveValue(row_to_select, "selected")
    #   }
    # })

    # Loop through sources to create server logic for DE and enrichment tables
    observeEvent(contrast(), {
      lapply(sources, function(source) {
        data_function <- reactive({ # this is not a function!! or, it is not a function in a reactive!!
          tryCatch(
            {
              if (is_filtering()) {
                req(row_to_select$selected)
                data <- movida_data$get_related_features(row_to_select$selected, source)
              } else {
                data <- movida_data$get_features_list(source)
              }

              data <- data.frame(features = data, stringsAsFactors = FALSE)
              rownames(data) <- data[, 1]


              if (!is.null(input[[paste0("feature_", source, "_search_box")]]) &&
                input[[paste0("feature_", source, "_search_box")]] != "") {
                search_term <- input[[paste0("feature_", source, "_search_box")]]
                data <- data[grep(search_term, rownames(data), ignore.case = TRUE), , drop = FALSE]
              }

              if (is.null(data) || length(data) == 0) {
                return(NULL)
              }

              return(data)
            },
            error = function(e) {
              warning(paste("Error in main_table_function:", e$message))
              return(NULL)
            }
          )
        })

        mod_table_server(
          id = paste0("feature_", source),
          main_table_function = data_function,
          selected_row = row_to_select,
          table_type = "minimal"
        )

        mod_table_server(
          id = paste0("de_", source),
          main_table_function = reactive({
            return(movida_data$get_dea(source, contrast(), FDRpvalue = NULL, FDRadj = NULL))
          }),
          selected_row = row_to_select
        )

        mod_table_server(
          id = paste0("enrich_", source),
          main_table_function = reactive({
            movida_data$get_fea(source, contrast())
          }),
          selected_row = row_to_select
        )
      })
    })
  })
}
