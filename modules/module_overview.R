library(shiny)
library(gargoyle)

source("modules/module_table_container.R")

# UI function for the module
mod_overview_ui <- function(id) {
  ns <- NS(id) # Namespace for the module

  # idea overview nuova.
  # le tabelle sono indipendent (sempre if)
  # tree pagine
  # tre div in tab organized \3 e ampio
  # dentro i div a search box (selecting), a list of ids, botton a button to is_filtering (reactive list if is_filtering)
  # seconda pagina -> tre tabelle, (o 2 o 1...)
  #
  # importante -> le tabelle mostrano solo i risultati disponibili, se no risultato, niente viene mostrato.
  # scelta contrast solo in de e enrichment
  # un altra cosa importante, devo poter mostrare piu plot modul at the same time

  # Sources for the module
  sources <- c(
    transcriptomics = "Transcriptomics",
    proteomics = "Proteomics",
    metabolomics = "Metabolomics"
  )

  contrast_selector_ui <- function(ns) {
    div(
      style = "display: flex; align-items: center; padding: 0 10px;",
      div(style = "margin-right: 10px;", "Select Contrast:"),
      div(style = "flex-grow: 1;", selectInput(ns("contrast_selector"), NULL, choices = NULL, width = "100%"))
    )
  }

  navset_pill(
    nav_panel(
      title = "Overview",
      layout_column_wrap(
        width = 1 / 3,
        lapply(names(sources), function(source) {
          div(
            h4(sources[[source]]),
            textInput(ns(paste0("feature_", source, "_search_box")), label = "Search:", placeholder = "Type to search..."),
            mod_table_ui(ns(paste0("feature_", source)))
          )
        })
      ),
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

    is_filtering <- reactiveVal(TRUE)

    # Update the contrast in the row_to_select reactive value
    observeEvent(contrast(), {
      row_to_select$contrast <- contrast()
      row_to_select$selected <- NULL
      row_to_select$source <- NULL
    })

    output$is_filtering_button <- renderUI({
      req(row_to_select$selected)
      if (is_filtering()) {
        actionButton(ns("is_filtering_button"), label = paste("Filtering based on", row_to_select$selected))
      } else {
        actionButton(ns("is_filtering_button"), label = paste("Click to filter based on", row_to_select$selected))
      }
    })

    # Observe the is_filtering button click to toggle is_filtering
    observeEvent(input$is_filtering_button, {
      is_filtering(!is_filtering())
    })

    sources <- c("transcriptomics", "proteomics", "metabolomics")

    observeEvent(is_filtering(), {
      if (is_filtering()) {
        freezeReactiveValue(row_to_select, "selected")
      }
    })

    # Loop through sources to create server logic for DE and enrichment tables
    lapply(sources, function(source) {
      mod_table_server(
        id = paste0("feature_", source),
        main_table_function = reactive({
          if (is_filtering()) {
            req(row_to_select$selected)
            data <- movida_data$get_related_features(row_to_select$selected, source)
          } else {
            data <- movida_data$get_features_list(source)
          }

          if (!is.null(input[[paste0("feature_", source, "_search_box")]]) &&
            input[[paste0("feature_", source, "_search_box")]] != "") {
            search_term <- input[[paste0("feature_", source, "_search_box")]]
            data <- data[grep(search_term, rownames(data), ignore.case = TRUE), , drop = FALSE]
          }

          if (is.null(data) || nrow(data) == 0) {
            return(NULL)
          }

          return(data)
        }),
        row_to_select = row_to_select
      )

      mod_table_server(
        id = paste0("de_", source),
        main_table_function = reactive({
          movida_data$get_de(source, contrast(), FDRpvalue = NULL, FDRadj = NULL)
        }),
        row_to_select = row_to_select
      )

      mod_table_server(
        id = paste0("enrich_", source),
        main_table_function = reactive({
          movida_data$get_enrichment(source, contrast())
        }),
        row_to_select = row_to_select
      )
    })
  })
}
