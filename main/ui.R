library(shiny)
library(bslib)
library(gargoyle)

source("modules/_module_dashboard.R")
source("modules/_module_overview.R")
source("modules/_module_sidebar_overview.R")

ui <- page_fillable(
  title = "MOVIDApp",
  theme = bs_theme(version = 5, brand = "_brand.yml"),

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
        sidebar = sidebar(uiOutput("sidebar_ui"), width = 350),  # dynamically rendered sidebar
        mod_overview_ui("overview")
      )
    )
  )
)