

server <- function(input, output, session, movida_data) {

  dashboard_elements <- reactiveValues(
    elements = list() # Initialize an empty list to store dashboard elements
  )

  selected_row_source <- reactiveValues(
    selected = NULL,
    source = NULL,
    contrast = NULL
  )

  # Reactive data frame to store and update data dynamically
  bookmarkedfeatures <- reactiveVal(data.frame())

  # UI
  # UI: Keep the button always present with a neutral initial state
  output$selectedfeature <- renderUI({
    req(selected_row_source$selected)

    feature <- isolate(safe_value(selected_row_source$selected))
    source <- isolate(safe_value(selected_row_source$source))
    contrast <- isolate(safe_value(selected_row_source$contrast))

    current_bookmarks <- isolate(bookmarkedfeatures())
    if (nrow(current_bookmarks) == 0) {
      current_bookmarks <- data.frame(
        feature = character(),
        source = character(),
        contrast = character(),
        stringsAsFactors = FALSE
      )
    }

    ui_info <- bookmark_ui_info(current_bookmarks, feature, source, contrast)

    div(
      style = "background-color: rgba(0, 128, 0, 0.8); height: 40px; display: flex; align-items: center; justify-content: space-between; padding: 20px 20px; color: white;",
      span(paste("Selected feature:", selected_row_source$selected, "from", selected_row_source$source)),
      div(
        style = "margin-left: auto;",
        actionButton(
          "bookmark_toggle",
          ui_info$label,
          style = sprintf("background-color: white; color: %s; border: none; padding: 5px 10px; height: 30px;", ui_info$color)
        )
      ),
      hr()
    )
  })


  observeEvent(input$bookmark_toggle,
    {
      feature <- isolate(safe_value(selected_row_source$selected))
      source <- isolate(safe_value(selected_row_source$source))
      contrast <- isolate(safe_value(selected_row_source$contrast))

      current_bookmarks <- isolate(bookmarkedfeatures())

      if (is_bookmarked(current_bookmarks, feature, source, contrast)) {
        updated_bookmarks <- remove_bookmark(current_bookmarks, feature, source, contrast)
      } else {
        updated_bookmarks <- add_bookmark(current_bookmarks, feature, source, contrast)
      }

      bookmarkedfeatures(updated_bookmarks)

      ui_info <- bookmark_ui_info(updated_bookmarks, feature, source, contrast)
      updateActionButton(session, "bookmark_toggle", label = ui_info$label)
      shinyjs::runjs(sprintf("$('#bookmark_toggle').css({'background-color':'white','color':'%s'});", ui_info$color))
    },
    ignoreInit = TRUE
  )




  # Observe and log the status of the bookmarkedfeatures variable
  observe({
    bookmarks <- bookmarkedfeatures()
    if (nrow(bookmarks) > 0) {
      print("Current bookmarks:")
      print(bookmarks)
    } else {
      print("No bookmarks present.")
    }
  })

  # Call module servers
  mod_overview_server("overview", dashboard_elements, selected_row_source, movida_data)
  mod_sidebar_overview_server("sidebar_overview", dashboard_elements, selected_row_source, movida_data)
  mod_dashboard_server("dashboard", dashboard_elements, movida_data) # Pass NULL or a list of elements as needed

  # For debugging purposes
  observe({
    # This observer will trigger whenever dashboard_elements$elements changes
    elements <- dashboard_elements$elements
    if (length(elements) > 0) {
      print("Current dashboard elements:")
      print(elements)
    } else {
      print("Dashboard elements list is empty.")
    }
  })

  sidebar_content <- reactiveVal(FALSE) # start with sidebar info hidden (set to FALSE)

  observeEvent(selected_row_source$contrast, {
    sidebar_content(FALSE) # Hide sidebar content when contrast changes (resetting the selection)
  })

  observeEvent(selected_row_source$selected, {
    sidebar_content(TRUE) # Hide sidebar content when contrast changes (resetting the selection)
  })

  # Render the sidebar UI based on the sidebar_content reactive value
  observeEvent(sidebar_content(), {
    output$sidebar_ui <- renderUI({
      if (sidebar_content()) {
        mod_sidebar_overview_ui("sidebar_overview")
      } else {
        # If sidebar is not visible, return an empty tagList
        tagList(h4("Select a row to show information"))
      }
    })
  })

  mod_plotting_server("plotting", dashboard_elements, bookmarkedfeatures, movida_data)
}