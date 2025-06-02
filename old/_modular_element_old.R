BaseElement <- R6Class("BaseElement",
  public = list(
    content = NULL,
    initialize = function(content) {
      self$content <- content
    },
    ui = function(id) {
      ns <- NS(id)

      btns <- list()
      controls <- list()

      if (!is.null(self$content$add_to_dashboard())) {
        btns <- append(btns, actionButton(ns("btn_dashboard"), "Add to Dashboard"))
      }
      if (!is.null(self$content$export_data())) {
        btns <- append(btns, actionButton(ns("btn_export_data"), "Export Data"))
      }
      if (!is.null(self$content$export_plot())) {
        btns <- append(btns, actionButton(ns("btn_export_plot"), "Export Plot"))
      }

      # If switch_visual returns a list of options, show selector
      visual_options <- self$content$switch_visual(NULL)
      if (!is.null(visual_options) && length(visual_options) > 0) {
        controls <- append(controls, selectInput(ns("select_visual"), "Select Visual", choices = visual_options))
      }

      tagList(
        div(class = "base-buttons", btns),
        div(class = "base-controls", controls),
        uiOutput(ns("main_output"))
      )
    },
    server = function(id, dashboard_list) {
      moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # For download content
        exported_data <- reactiveVal(NULL)
        exported_plot <- reactiveVal(NULL)

        observeEvent(input$btn_dashboard, {
          result <- self$content$add_to_dashboard()
          if (!is.null(result)) {
            key <- paste0("elem_", as.character(Sys.time()))
            dashboard_list[[key]] <- result
          }
        })

        observeEvent(input$btn_export_data, {
          data <- self$content$export_data()
          exported_data(data)

          showModal(modalDialog(
            title = "Download Data",
            downloadButton(ns("download_data"), "Download CSV"),
            easyClose = TRUE
          ))

          output$download_data <- downloadHandler(
            filename = function() paste0("data-", Sys.Date(), ".csv"),
            content = function(file) {
              write.csv(exported_data(), file, row.names = FALSE)
            }
          )
        })

        observeEvent(input$btn_export_plot, {
          plot <- self$content$export_plot()
          exported_plot(plot)

          showModal(modalDialog(
            title = "Download Plot",
            downloadButton(ns("download_plot"), "Download PNG"),
            easyClose = TRUE
          ))

          output$download_plot <- downloadHandler(
            filename = function() paste0("plot-", Sys.Date(), ".png"),
            content = function(file) {
              png(file, width = 800, height = 600)
              print(exported_plot())
              dev.off()
            }
          )
        })

        observeEvent(input$select_visual,
          {
            selected <- input$select_visual
            result_plot <- self$content$switch_visual(selected)

            if (!is.null(result_plot)) {
              output$main_output <- renderUI({
                plotOutput(ns("visual_plot"))
              })
              output$visual_plot <- renderPlot({
                result_plot
              })
            }
          },
          ignoreInit = TRUE
        )
      })
    }
  )
)
