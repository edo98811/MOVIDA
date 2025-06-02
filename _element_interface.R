library(R6)

MainContentInterface <- R6Class("MainContentInterface",
  public = list(

    #' Add visual to the dashboard
    add_to_dashboard = function() {
      NULL
    },
    
    #' Export data as csv
    export_data = function() {
      NULL
    },
    
    #' Export plot as png
    export_plot = function() {
      NULL
    },
    
    #' Switch internal visual state (e.g., change plot type)
    switch_visual = function(view) {
      NULL
    }
  )
)
