mod_sidebar_overview_ui <- function(id) {
  ns <- NS(id)
  # here i want a navset pill with the three type of omics, i want to see as many plots as there are related features for each type
  # i want a selector for just related contrast and for the column of the metadata to use to plot (maybe son)


  sources <- c(
    transcriptomics = "Transcriptomics",
    proteomics = "Proteomics",
    metabolomics = "Metabolomics"
  )

  div(
    uiOutput(ns("sidepanel_title")),
    do.call(
      navset_pill,
      lapply(names(sources), function(source) {
        nav_panel(
          title = sources[[source]],
          div(
            uiOutput(ns(paste0("sidebar_ui_plots_", source)))
          )
        )
      })
    )
  )
}

mod_sidebar_overview_server <- function(id, dashboard_elements, selected_row_source, movida_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    sources <- movida_data$get_sources()

    # sidebar with the heatmaps
    # is it maybe not necessary? can i keep it in the ui?
    observeEvent(selected_row_source$selected, {
      if (suppressWarnings(check_goterm(selected_row_source$selected))) {
        #
        titles <- c(
          transcriptomics = "Pathway Genes Expression",
          proteomics = "Protein Abundance",
          metabolomics = "Metabolite Abundance"
        )
        # sidebar with the heatmaps
        lapply(sources, function(source) {
          output[[paste0("sidebar_ui_plots_", source)]] <- renderUI({
            req(selected_row_source$selected)
            pathway <- selected_row_source$selected

            # return UI elements for filtering
            # div(
            #   class = "filter-controls mb-3",
            #   selectInput(
            #     ns(paste0("group_column_", source)),
            #     "Group by:",
            #     choices = NULL,
            #     selected = "condition"
            #   ),
            #   uiOutput(ns(paste0("selected_contrasts_only_", source))),
            #   class = "btn-group", role = "group",
            #   tags$button(id = ns(paste0("bookmark_all", pathway)), type = "button", class = "btn btn-primary", "Bookmark")
            # ),
            mod_plot_ui(
              id = ns(paste0(source, "_heatmap_", gsub(":", "", pathway))),
              title = titles[[source]],
              show_export_data_btn = TRUE,
              show_export_plot_btn = TRUE
            )
          })
        })

        # SERVER LOGIC -------
        lapply(sources, function(source) {
          req(selected_row_source)

          pathway <- selected_row_source$selected
          # source <- selected_row_source$source
          contrast <- selected_row_source$contrast

          # Plot id for the heatmap
          plot_id <- reactive({
            paste0("source_", source, ";feature_", pathway)
          })

          plot_function <- function(export_data = FALSE) {
            features <- movida_data$get_pathwayfeatures(pathway, contrast, source)
            se <- movida_data$get_values_all(source, return_se = TRUE)

            plot_heatmap_movida(
              se,
              features,
              row_column = "SYMBOL",
              export_data = export_data
            )
          }

          # Set up plot server logic for this feature
          mod_plot_server(
            id = paste0(source, "_heatmap_", gsub(":", "", pathway)),
            plot_id = plot_id(),
            main_plotting_function = plot_function,
            dashboard_elements = dashboard_elements
          )
        })
      } else if (
        check_is_valid_feature(selected_row_source$selected)) {
        titles <- c(
          transcriptomics = "Gene Expression",
          proteomics = "Protein Abundance",
          metabolomics = "Metabolite Abundance"
        )
        # sidebar with geneplots

        # Module: Sidebar Overview
        #
        # This module creates a tabbed interface (pill) with 3 tabs. It's responsible for:
        # 1. Fetching related features based on the user's currently selected feature
        # 2. Generating a plot for each related feature
        # 3. Providing interactive elements including:
        #    - A button to bookmark features
        #    - A button to select features
        #
        # The module handles returning either the selected feature or all related features
        # depending on the context.

        # sidebar with the gene plots
        # RENDER UI -------
        # 1. For each omics source (transcriptomics, proteomics, metabolomics):
        #    - Add buttons to filter by contrast (maybe another select contrast specific to the navbar overview) and to select grouping columns
        # 2. When a feature is selected:
        #    - Retrieve related features using get_related_features()
        #    - Generate a plot for each related feature
        #    - Add selection and bookmark buttons for each plot
        # 3. Update UI reactively when selection changes

        # It is legal to put renderUI inside lapply: https://groups.google.com/g/shiny-discuss/c/DFcPtnpfETs (joe cheng)
        # Iterate through each source and create a UI output for each
        output$sidepanel_title <- renderUI({
          req(selected_row_source$selected)
          return(tags$h4(paste0("features related to ", selected_row_source$selected)))
        })

        lapply(sources, function(source) {
          output[[paste0("sidebar_ui_plots_", source)]] <- renderUI({
            req(selected_row_source$selected)
            relatedfeatures <- movida_data$get_related_features(selected_row_source$selected, source)

            # return UI elements for filtering
            list(
              div(
                class = "filter-controls mb-3",
                br(),
                selectInput(
                  ns(paste0("group_column_", source)),
                  "Group by:",
                  choices = c("group"),
                  selected = "group"
                ),
                # Source selector
                selectInput(
                  ns(paste0("subset_group_", source)),
                  "Select Subset",
                  choices = NULL,
                  selected = NULL,
                  multiple = TRUE
                )
              ),
              # Generate a plot for each related feature
              lapply(relatedfeatures, function(feature) {
                div(
                  hr(),
                  mod_plot_ui(
                    id = ns(paste0(source, "_plot_", feature)),
                    title = titles[[source]],
                    show_export_data_btn = TRUE,
                    show_export_plot_btn = TRUE
                  ),
                  div(
                    class = "btn-group w-100 py-2", role = "group",
                    tags$button(id = ns(paste0("select_", feature)), type = "button", class = "btn btn-outline-secondary btn-sm", "Select"),
                    tags$button(id = ns(paste0("bookmark_", feature)), type = "button", class = "btn btn-outline-secondary btn-sm", "Bookmark")
                  )
                )
              })
            )
          })
        })

        # SERVER LOGIC -------
        # Module for handling the sidebar overview functionality in MOVIDA project.
        # This code orchestrates the interaction with various data sources by:
        # 1. Iterating through each data source
        # 3. Processing all features related to the current source by:
        #    - Creating bookmark logic (if applicable)
        #    - Implementing selection logic
        #    - Generating plotting functions
        #    - Creating unique plot IDs
        #    - Setting up server-side module functionality


        lapply(sources, function(source) {
          req(selected_row_source$selected)
          relatedfeatures <- movida_data$get_related_features(selected_row_source$selected, source)

          # observe event for select filtering column (does this need to be rectie?)
          group_by <- reactive({
            input[[paste0("group_column_", source)]]
          })

          observeEvent(input[[paste0("group_column_", source)]], {
            input_id <- paste0("group_column_", source)

            # Only update if the input has already been rendered
            if (!is.null(input[[input_id]])) {
              updateSelectInput(
                session,
                input_id,
                choices = colnames(movida_data$get_metadata(source)),
                selected = group_by()
              )
            }
          })

          # Observe the group column selection and update the subset group choices
          observe({
            selected_group_col <- input[[paste0("group_column_", source)]]
            req(selected_group_col)

            updateSelectInput(
              session,
              paste0("subset_group_", source),
              choices = movida_data$get_metadata(source)[[selected_group_col]],
              selected = "group"
            )
          })

          # For each related feature, set up the select and bookmark button observers
          # Then create the plotting function, call the plot_module server
          lapply(relatedfeatures, function(feature) {
            selected_metadata_subset <- reactive({
              input[[paste0("subset_group_", source)]]
            })
            selected_metadata_column <- reactive({
              input[[paste0("group_column_", source)]]
            })

            # Create observer for the select button for this feature
            observeEvent(input[[paste0("select_", feature)]], {
              selected_row_source$selected <- feature
              selected_row_source$source <- source
            })

            # Create observer for the bookmark button for this feature
            observeEvent(input[[paste0("bookmark_", feature)]], {
              # Bookmark functionality to be implemented
            })

            plot_id <- reactive({
              paste0(
                "source_", source, ";selected_", feature,
                ";selected_metadata_subset", selected_metadata_subset(),
                ";selected_metadata_column", selected_metadata_column()
              )
            })

            # Generate the plot function for this feature
            # I don't need it reactive because it is created new each time
            plot_function <- reactive({
              function(export_data = FALSE) {
                plot_expression_movida(
                  feature,
                  movida_data$get_values_subset_metadata(selected_metadata_subset(), source, column = selected_metadata_column(), return_se = TRUE),
                  export_data = export_data,
                  data_type = "source",
                  group_var = group_by()
                )
              }
            })

            # Set up plot server logic for this feature
            mod_plot_server(
              id = paste0(source, "_plot_", feature),
              plot_id = plot_id(),
              main_plotting_function = plot_function(),
              dashboard_elements = dashboard_elements
            )
          })
        })
      }
    })
  })
}
