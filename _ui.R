library(shiny)
library(bslib)
library(gargoyle)

source("1_module_dashboard.R")
source("1_module_overview.R")
source("1_module_sidebar_overview.R")

ui <- page_fillable(
  title = "MOVIDApp",
  theme = bs_theme(version = 5, brand = "_brand.yml"),
  head = tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "sidebar-toggle.css"),
    tags$script(src = "sidebar-toggle.js")
  ),

  page_navbar(
    title = "MOVIDA",
    navbar_options = navbar_options(
      underline = TRUE
    ),

    nav_panel(
      title = "Dashboard",
      mod_dashboard_ui("dashboard")
    ),

    nav_panel(
      title = "Overview",
      layout_sidebar(
        sidebar = sidebar(uiOutput("sidebar_ui")),  # dynamically rendered sidebar
        mod_overview_ui("overview")
      )
    )
  )
)