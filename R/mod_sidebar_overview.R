mod_sidebar_overview_ui <- function(id) {
  ns <- NS(id)
  # here i want a navset pill with the three type of omics, i want to see as many plots as there are related features for each type
  # i want a selector for just related contrast and for the column of the metadata to use to plot (maybe son)

  div(
    uiOutput(ns("sidepanel_title")),
    uiOutput(ns("sidebar_ui_plots")),
  )
}
mod_sidebar_overview_server <- function(id, dashboard_elements, selected_row_source, movida_data, max_plots = 20) {
  moduleServer(id, function(input, output, session) {
     ns <- session$ns
    sources <- movida_data$get_sources()

    # Store current features per source
    current_features <- reactiveValues(
      list = setNames(vector("list", length(sources)), names(sources))
    )

    # Update current features when selected row changes
    observeEvent(selected_row_source$selected, {
      req(selected_row_source$selected)
      cat("Selected feature changed to:", selected_row_source$selected, "\n")

      for (source in names(sources)) {
        cat("Getting related features for source:", source, "\n")
        related_features <- movida_data$get_related_features(
          selected_row_source$selected, selected_row_source$source,
          target = source
        )
        cat("Found", length(related_features %||% c()), "related features\n")
        current_features$list[[source]] <- related_features
      }
    })

    # Sidebar navigation
    output$sidebar_ui_plots <- renderUI({
      cat("Rendering sidebar UI, selected:", selected_row_source$selected, "\n")

      do.call(
        navset_pill,
        lapply(names(sources), function(source) {
          nav_panel(
            title = sources[[source]],
            tagList(
              br(),
              tags$h4(selected_row_source$selected),
              uiOutput(ns(paste0("sidebar_ui_plots_", source)))
            )
          )
        })
      )
    })

    # Per-source UI
    lapply(names(sources), function(source) {
      output[[paste0("sidebar_ui_plots_", source)]] <- renderUI({
        features <- current_features$list[[source]]

        # Debug: Print features info to console
        cat("Rendering UI for source:", source, "\n")
        cat("Features found:", length(features %||% c()), "\n")
        if (!is.null(features)) {
          cat("Feature names:", paste(features, collapse = ", "), "\n")
        }

        tagList(
          # Filters
          div(
            class = "filter-controls mb-3",
            selectInput(
              ns(paste0("group_column_", source)),
              "Group by:",
              choices = movida_data$get_metadata_columns(source),
              selected = movida_data$get_metadata_columns(source)[[1]]
            ),
            selectInput(
              ns(paste0("subset_group_", source)),
              "Select Subset",
              choices = NULL,
              selected = NULL,
              multiple = TRUE
            )
          ),
          # Only render plots for existing features
          if (!is.null(features) && length(features) > 0) {
            lapply(seq_len(min(length(features), max_plots)), function(i) {
              if (i <= length(features) && !is.null(features[[i]])) {
                feature <- features[[i]]
                div(
                  hr(),
                  mod_plot_ui(
                    id = ns(paste0(source, "_plot_", i)),
                    title = feature,
                    show_export_data_btn = TRUE,
                    show_export_plot_btn = TRUE
                  ),
                  div(
                    class = "btn-group w-100 py-2",
                    tags$button(
                      id = ns(paste0("select_", source, "_", i)),
                      type = "button",
                      class = "btn btn-outline-secondary btn-sm",
                      "Select"
                    )
                  )
                )
              }
            })
          } else {
            # Show placeholder when no features available
            div(
              class = "text-center text-muted p-4",
              "No related features found for this source"
            )
          }
        )
      })
    })

    # Initialize all plot modules
    lapply(names(sources), function(source) {
      # Update subset group choices
      observeEvent(        
        {
          input[[paste0("group_column_", source)]]
          selected_row_source$selected
        }, {
        selected_group_col <- input[[paste0("group_column_", source)]]
        req(selected_group_col)
        updateSelectInput(
          session,
          paste0("subset_group_", source),
          choices = movida_data$get_metadata(source)[[selected_group_col]],
            selected = head(movida_data$get_metadata(source)[[selected_group_col]], 5)
        )
      })

      lapply(seq_len(max_plots), function(i) {
        plot_id <- reactive({
          features <- current_features$list[[source]]
          # Check if features is NULL, if features is shorter
          if (is.null(features) || length(features) < i || is.null(features[[i]])) {
            return(NULL)
          }
          paste0("source_", source, ";feature_", features[[i]])
        })

        # The function with the ractive values inside is passed to the plotting
        main_plot_function <- function(export_data = FALSE) {
          features <- current_features$list[[source]]

          if (is.null(features) || length(features) < i || is.null(features[[i]])) {
            return(NULL)
          }

          feature <- features[[i]]
          group <- input[[paste0("subset_group_", source)]]
          column <- input[[paste0("group_column_", source)]]

          return(plot_expression_movida(
            feature,
            movida_data$get_values_subset_metadata(
              group,
              source,
              column = column,
              return_se = TRUE
            ),
            export_data = export_data,
            data_type = "source",
            group_var = column
          ))
        }

        # start the feature mod plot server
        mod_plot_server(
          id = paste0(source, "_plot_", i),
          plot_id = plot_id,
          main_plotting_function = main_plot_function,
          dashboard_elements = dashboard_elements
        )

        # Select
        observeEvent(input[[paste0("select_", source, "_", i)]], {
          features <- current_features$list[[source]]
          if (!is.null(features) && length(features) >= i && !is.null(features[[i]])) {
            selected_row_source$selected <- features[[i]]
            selected_row_source$source <- source
          }
        })
      })
    })
  })
}
